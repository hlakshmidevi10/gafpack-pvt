use std::cmp::min;
use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{prelude::*, BufReader, BufWriter};
use std::path::Path;
use std::io::Write;
use std::time::Instant;
use bytemuck::{Pod, Zeroable};
use rayon::prelude::*;

/// On-disk record in `_path_pos_v2.bin` (written by find_mems v2). Little-endian, 16 bytes.
///
/// See `mem-projection/pangenome-pipeline/PLAN_find_mems_binary_io_v2.md`.
/// node_id and offset are NOT stored: gafpack derives them by walking the
/// path's cum_bp prefix-sum and locating the step containing path_bp.
/// Records are sorted by path_bp within each seq_id bucket so the walker
/// can advance a single monotonic cursor (linear merge).
///
/// is_rev(graph_pos) from the BWT hit is NOT stored either. The GAF strand
/// is the XOR of the bucket's path-strand (seq_id & 1) and the GFA step's
/// +/- orientation; see traverse_nodes(). An earlier v2 draft stored a
/// "sanity bit" in match_len bit 31 and asserted it equalled the bucket
/// parity -- that conflated two independent strand concepts and silently
/// dropped valid records. Don't reintroduce it.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Record {
    path_bp: u32,
    match_len: u32,
    read_st: u32,
    read_id: u32,
}

// Compile-time guarantee that the on-disk size stays 16 bytes.
const _: () = assert!(std::mem::size_of::<Record>() == 16);

/// Iterates through each line in a file, applying the provided callback function
///
/// # Arguments
/// * `filename` - Path to the file to read
/// * `callback` - Function to call for each line
fn for_each_line_in_file(filename: &str, mut callback: impl FnMut(&str)) {
    let file = File::open(filename).unwrap();
    let (reader, _compression) = niffler::get_reader(Box::new(file)).unwrap();
    let buf_reader = BufReader::new(reader);
    for line in buf_reader.lines() {
        callback(&line.unwrap());
    }
}

/// Process each step in a GAF alignment line, calculating coverage for graph nodes
///
/// # Arguments
/// * `line` - A GAF format alignment line
/// * `callback` - Function called for each node with (node_id, coverage_length)
/// * `get_node_len` - Function to get the length of a node by its ID
///
/// # Details
/// Parses GAF alignment lines to extract node coverage information:
/// - Handles both forward (>) and reverse (<) node traversals
/// - Adjusts coverage for partial node alignments at path ends
/// - Accumulates coverage across multi-node paths
fn for_each_step(
    line: &str,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
) {
    //eprintln!("{}", line);
    let walk = line.split('\t').nth(5).unwrap();
    if walk != "*" {
        //eprintln!("oheunotoeunthoue");
        let target_start = line.split('\t').nth(7).unwrap().parse::<usize>().unwrap();
        let target_end = line.split('\t').nth(8).unwrap().parse::<usize>().unwrap();
        let target_len = target_end - target_start;
        //eprintln!("target_len = {}", target_len);
        let fields = line
            .split('\t')
            .nth(5)
            .unwrap()
            .split(|c| c == '<' || c == '>')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<usize>().unwrap())
            .enumerate()
            .collect::<Vec<(usize, usize)>>();
        let mut seen: usize = 0;
        let fields_len = fields.as_slice().len();
        //eprintln!("fields len = {}", fields_len);
        for (i, j) in fields {
            let mut len = get_node_len(j);
            //eprintln!("node {} len = {}", j, len);
            if i == 0 {
                //eprintln!("on first step {} {} {}", len, target_start, seen);
                assert!(len >= target_start);
                len -= target_start;
            }
            if i == fields_len - 1 {
                //eprintln!("on last step {} {} {}", len, target_end, seen);
                assert!(target_len >= seen);
                len = target_len - seen;
            }
            if i == fields_len {
                assert!(false);
            }
            //eprintln!("node {} adj len = {}", j, len);
            seen += len;
            callback(j, len);
        }
        //eprintln!("seen = {}", seen);
    }
    //eprintln!("at end");
}

/// Create a reader that handles compressed files
fn create_reader(path: &Path) -> std::io::Result<Box<dyn BufRead>> {
    let file = File::open(path)?;

    // 4 MiB rather than BufReader's 8 KiB default. The GFA is scanned end to
    // end twice (S-lines in parse_gfa, P-lines in the walker) and reaches
    // 26.4 GB on HPRCv2.1 MC chr1 -- at 8 KiB that is ~3.3M read syscalls per
    // pass. Sequential scan, so a large buffer costs nothing but the pages.
    const GFA_BUF: usize = 4 << 20;
    if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::with_capacity(GFA_BUF, decoder);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(GFA_BUF, file);
        Ok(Box::new(buf_reader))
    }
}


/// Strip a trailing CR/LF from a raw line.
#[inline]
fn trim_eol(b: &[u8]) -> &[u8] {
    let mut e = b.len();
    while e > 0 && (b[e - 1] == b'\n' || b[e - 1] == b'\r') {
        e -= 1;
    }
    &b[..e]
}

/// ASCII-decimal parse straight off bytes, so the GFA scan never has to
/// materialise a &str. Returns None on empty input or any non-digit.
#[inline]
fn parse_usize_bytes(b: &[u8]) -> Option<usize> {
    if b.is_empty() {
        return None;
    }
    let mut v: usize = 0;
    for &c in b {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as usize)?;
    }
    Some(v)
}

/// Parse GFA file and extract segment information
/// Returns (segment_lengths, min_id) where segment_lengths[id - min_id] gives the length
fn parse_gfa(gfa_path: &str) -> std::io::Result<(Vec<usize>, usize)> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    // Read raw bytes rather than read_line-into-String: the GFA reaches
    // 26.4 GB on HPRCv2.1 MC chr1 and read_line UTF-8-validates every byte of
    // it. Nothing here needs a &str -- ids and sequence lengths come straight
    // off the bytes.
    let mut buf: Vec<u8> = Vec::with_capacity(1 << 16);
    // min_id isn't known until the full scan completes (segment lines can
    // arrive in any id order), so the dense array can't be sized/indexed
    // during the scan itself. Collect into a plain append-only Vec instead
    // of a HashMap -- no hashing at all (~4.7M inserts on HPRC chr6) -- then
    // do a second short pass to fill the dense array once min_id is known.
    let mut segments: Vec<(usize, usize)> = Vec::new();
    let mut min_id = usize::MAX;
    let mut max_id = 0;

    loop {
        buf.clear();
        let bytes_read = reader.read_until(b'\n', &mut buf)?;
        if bytes_read == 0 {
            break;
        }

        let line_b = trim_eol(&buf);

        // Only process segment lines
        if line_b.first() != Some(&b'S') {
            continue;
        }

        // Parse segment line format: S<tab>id<tab>sequence
        let mut fields = line_b.split(|&c| c == b'\t');
        let Some((id_b, seq)) = fields.next().and_then(|_type| {
            let id_b = fields.next()?;
            let seq = fields.next()?;
            Some((id_b, seq))
        }) else {
            continue;
        };

        // Parse segment ID
        let id = parse_usize_bytes(id_b).expect("non-numeric GFA segment id");
        min_id = min_id.min(id);
        max_id = max_id.max(id);
        segments.push((id, seq.len()));
    }

    // Create a dense vector for O(1) access
    let num_segments = max_id - min_id + 1;
    let mut segment_lengths = vec![0; num_segments];

    for (id, len) in segments {
        segment_lengths[id - min_id] = len;
    }

    Ok((segment_lengths, min_id))
}


/// Advance a monotonic step cursor so that `cum_bp[*i] <= path_bp < cum_bp[*i+1]`.
///
/// Pure helper, extracted from `process_path_matches` so the derivation can be
/// unit-tested in isolation. Called once per record in the v2 walker. Because
/// records are sorted by `path_bp` ascending within a seq_id bucket, the cursor
/// only advances forward across the whole bucket -- amortized O(1) per call.
///
/// Preconditions (debug-assert only; production callers guarantee them):
/// - cum_bp.len() >= 2 (path has >= 1 step)
/// - cum_bp[0] == 0 and cum_bp is non-decreasing
/// - path_bp < *cum_bp.last() (caller filters path_bp >= path_total_bp)
/// - *i < cum_bp.len() - 1
/// - cum_bp[*i] <= path_bp (cursor never goes backwards)
#[inline]
fn advance_step_cursor(
    path_bp: usize,
    cum_bp: &[usize],
    i: &mut usize,
    cum: &mut usize,
) {
    let last_step = cum_bp.len() - 1; // == n_steps; cum_bp has n_steps+1 entries
    while *i + 1 < last_step && path_bp >= cum_bp[*i + 1] {
        *i += 1;
        *cum = cum_bp[*i];
    }
}

fn traverse_nodes(
    steps: &[&str],
    node_ids: &[usize],
    callback: &mut impl FnMut(usize, usize),
    get_node_len: &mut impl FnMut(usize) -> usize,
    start_offset: usize,
    match_len: usize,
    is_reverse: bool,
    path_str: &mut String,
) -> usize
{
    if steps.is_empty() || match_len == 0 {
        return 0;
    }

    // GAF strand is the XOR of the node's GFA orientation (+/-) and the path
    // walk direction. is_reverse is constant for this whole call, so only 2
    // strand values are ever possible here -- precompute both once (stack
    // values) instead of re-deriving via a match on every step-visit.
    let (plus_strand, minus_strand) = if is_reverse { ('<', '>') } else { ('>', '<') };

    let mut total_path_length = 0;
    let mut remaining_len = match_len;

    // zip, not enumerate()+index: node_id per step precomputed once per path
    // by the caller (see process_path_matches) to avoid re-parsing digits on
    // every step-visit (called once per MEM record touching this step;
    // ~805M step-visits measured on HPRC chr6). zip's bounds tracking is
    // inherent to each iterator's own exhaustion, not an extra cross-slice
    // length check the way index-based access (steps[i], node_ids[i]) is.
    for (i, (step, &node_id)) in steps.iter().zip(node_ids.iter()).enumerate() {
        let (seg, orient) = step.split_at(step.len() - 1);
        let node_length = get_node_len(node_id);
        total_path_length += node_length;

        let strand = if orient == "+" { plus_strand } else { minus_strand };
        path_str.push(strand);
        path_str.push_str(seg);

        let coverage = if i == 0 {
            // First node: account for start_offset
            min(remaining_len, node_length - start_offset)
        } else {
            // Subsequent nodes: use full node length or remaining length
            min(remaining_len, node_length)
        };

        remaining_len -= coverage;
        callback(node_id, coverage);

        if remaining_len == 0 {
            break;
        }
    }

    total_path_length
}

/// Read path positions CSV file and return the data
fn read_path_pos_csv(path_pos_file: &str) -> std::io::Result<Vec<(usize, usize, u8, usize, usize, String)>> {
    let mut data = Vec::new();
    let file = File::open(path_pos_file)?;
    let reader = BufReader::new(file);

    for (line_num, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();

        // Skip header line if present
        if line_num == 0 && (line.contains("node_id") || line.contains("offset")) {
            continue;
        }

        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split_whitespace().collect();
        if fields.len() >= 6 {
            let node_id = fields[0].trim().parse::<usize>().unwrap();
            let offset = fields[1].trim().parse::<usize>().unwrap();
            let strand = fields[2].trim().parse::<u8>().unwrap();
            let mem_len = fields[3].trim().parse::<usize>().unwrap();
            let read_st = fields[4].trim().parse::<usize>().unwrap();
            let read_id = fields[5].trim().to_string();

            data.push((node_id, offset, strand, mem_len, read_st, read_id));
        }
    }

    Ok(data)
}

/// Read and process seq_id_starts file
/// Returns a HashMap mapping seq_id -> (cumulative_start, count)
///
/// File format: seq_id\tcumulative_length
/// Special case: -1 indicates end of processing
///
/// Example input:
/// 0    0
/// 1    3
/// 2    5
/// -1   10
///
/// Results in:
/// 0 -> (0, 3)
/// 1 -> (3, 2)
/// 2 -> (5, 5)
fn read_seq_id_starts_map(seq_id_starts_file: &str) -> std::io::Result<HashMap<String, (usize, usize)>> {
    let mut seq_map = HashMap::new();
    let file = File::open(seq_id_starts_file)?;
    let reader = BufReader::new(file);

    let mut lines: Vec<(String, usize)> = Vec::new();

    // Read all lines first
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();

        // eprintln!("Line: {}", line);

        if line.is_empty() {
            continue;
        }
        // eprintln!("Line bytes: {:?}", line.as_bytes());

        let fields: Vec<&str> = line.split(' ').collect();
        if fields.len() >= 2 {
            let seq_id = fields[0].trim().to_string();
            let cumulative_length = fields[1].trim().parse::<usize>().unwrap();
            lines.push((seq_id, cumulative_length));
        }
    }

    // Process the lines to calculate counts
    for i in 0..lines.len() {
        let (seq_id, cumulative_start) = &lines[i];

        // Skip the -1 terminator
        if seq_id == "-1" {
            break;
        }

        let count = if i + 1 < lines.len() {
            let (next_seq_id, next_cumulative) = &lines[i + 1];
            if next_seq_id == "-1" {
                // Last entry uses the -1 line's cumulative value as the total
                next_cumulative - cumulative_start
            } else {
                next_cumulative - cumulative_start
            }
        } else {
            // If there's no next line, count is 0 (shouldn't happen with proper format)
            0
        };

        eprintln!("Storing: seq_id: {}, val: {}, {}", seq_id, cumulative_start, count);

        seq_map.insert(seq_id.clone(), (*cumulative_start, count));
    }

    Ok(seq_map)
}

/// Reads a file containing cumulative start values, returning a vector
/// where vec[seq_id] = cumulative_start
///
/// # Arguments
/// * `filename` - Path to the file containing cumulative start values
///
/// # File Format
/// Each line contains a single cumulative start value.
/// The sequence ID is implied by the line number (0-indexed).
/// Example:
/// ```
/// 1000
/// 2500
/// 4000
/// ```
/// Results in: vec[0] = 1000, vec[1] = 2500, vec[2] = 4000
///
/// # Returns
/// A vector where the index is the sequence ID and value is the cumulative start position
fn read_cumulative_starts(filename: &str) -> std::io::Result<Vec<usize>> {
    let mut cumulative_starts = Vec::new();

    for_each_line_in_file(filename, |line: &str| {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return;
        }

        if let Ok(cumulative_start) = line.parse::<usize>() {
            cumulative_starts.push(cumulative_start);
        }
    });

    // eprintln!("Cumulative starts vector:");
    // for (index, value) in cumulative_starts.iter().enumerate() {
    //     eprintln!("  [{}]: {}", index, value);
    // }
    eprintln!("Size of cumulative starts: {}; Total no. of entries: {}", cumulative_starts.len(), cumulative_starts[cumulative_starts.len() - 1]);

    Ok(cumulative_starts)
}

/// v2 path-walker: linear merge of records (sorted by path_bp ascending)
/// against the path's steps (advanced via a monotonic cursor).
///
/// Invariant per iteration: cum_bp[i] <= r.path_bp < cum_bp[i+1].
/// Cursor `i` is fresh per call -- the forward and reverse orientations get
/// disjoint record slices (different seq_ids) and so disjoint walker state.
/// Total work across the slice: O(steps + records), no per-record search.
///
/// Replaces the v1 cum_bp.partition_point() lookup. Also removes the v1
/// mismatch_node/mismatch_off counters -- their underlying redundancy
/// (storing node_id/offset alongside path_bp) is gone in v2.
fn process_path_matches(
    steps: &[&str],
    seq_id: usize,
    name: &str,
    is_reverse: bool,
    path_starts: &[usize],
    records: &[Record],
    gaf_output: &mut Option<BufWriter<File>>,
    callback: &mut impl FnMut(usize, usize),
    get_node_len: &mut impl FnMut(usize) -> usize,
    // When Some, dedup records by (read_id, read_st, starting_node_id, offset)
    // packed into a single u64 -- see walk_gfa() for the bit layout and the
    // load-time field-width asserts. Each unique 4-tuple contributes coverage
    // and GAF output exactly once across the whole gafpack run; subsequent
    // records with the same key are skipped. See the --dedup-read-node CLI
    // flag for rationale.
    dedup_seen: &mut Option<HashSet<u64>>,
    dedup_skipped: &mut u64,
    // Minimum raw GFA node_id; subtracted from the starting node before
    // packing into the dedup key so the 25-bit field stores a min_id-relative
    // value (range 0..=max_id-min_id) rather than the raw id. Some graph
    // builds (e.g. cactus-mc) use sparse node-id ranges where raw max_id
    // overflows 25 bits even though (max_id - min_id) fits comfortably.
    min_node_id: usize,
    verbose: bool,
) -> std::io::Result<(usize, u64)> {

    let st = path_starts[seq_id];
    let end = path_starts[seq_id + 1];
    let num_matches = end - st;

    if num_matches == 0 {
        return Ok((0, 0));
    }

    if verbose {
        eprintln!("--------");
        eprintln!("Processing path: {}   seq_id: {}   st: {}   num_matches: {}", name, seq_id, st, num_matches);
    }

    // Pre-parse steps once: node_id + GAF strand char per step, plus cumulative
    // bp prefix sum. Both node_id and strand are invariant per (path, step,
    // is_reverse) -- traverse_nodes previously re-derived them (re-parsing the
    // digits and re-running the strand match) on every one of the ~805M
    // step-visits measured on HPRC chr6, once per MEM record touching that
    // step. Precomputing here amortizes the cost to once per path.
    // Pre-parse steps once: node_id per step + cumulative bp prefix sum.
    // (Strand is NOT precomputed into a parallel array -- an earlier attempt
    // at that regressed wall time on vesuvio despite removing work, most
    // likely from the added cross-array bounds-checked indexing outweighing
    // the eliminated match. is_reverse is constant for the whole call, so
    // traverse_nodes derives strand from 2 precomputed stack constants and
    // the orient byte it already has for free, instead of a stored array.)
    let mut step_node_ids: Vec<usize> = Vec::with_capacity(steps.len());
    let mut cum_bp: Vec<usize> = Vec::with_capacity(steps.len() + 1);
    cum_bp.push(0);
    for step in steps.iter() {
        let (seg, _orient) = step.split_at(step.len() - 1);
        let seg_id = seg.parse::<usize>().unwrap();
        step_node_ids.push(seg_id);
        cum_bp.push(cum_bp.last().unwrap() + get_node_len(seg_id));
    }
    let path_total_bp = *cum_bp.last().unwrap();
    let n_steps = step_node_ids.len();
    debug_assert!(n_steps > 0, "path with zero steps");

    let mut total_gaf_entries = 0;
    let mut path_str = String::with_capacity(64);

    // Sanity counter for out-of-range path_bp. Should be zero on healthy runs.
    let mut path_bp_out_of_range: u64 = 0;

    // Monotonic step cursor: i = current step index, cum = cum_bp[i].
    // Records are sorted by path_bp ascending within this bucket, so i
    // only ever advances forward. Amortized O(1) per record.
    let mut i: usize = 0;
    let mut cum: usize = 0;

    for idx in st..end {
        let r = &records[idx];
        let path_bp = r.path_bp as usize;
        let match_len = r.match_len as usize;
        let read_start = r.read_st as usize;
        let read_id = r.read_id as usize;

        // path_bp >= path_total_bp: hit lands past the end of the path.
        // Preserved from v1; same behavior, same threshold (>=, not >).
        if path_bp >= path_total_bp {
            path_bp_out_of_range += 1;
            if path_bp_out_of_range <= 5 {
                eprintln!("ERROR: path_bp {} >= path length {} on path {} (seq_id {}); skipping record idx {}",
                          path_bp, path_total_bp, name, seq_id, idx);
            }
            continue;
        }

        // Advance the cursor so cum_bp[i] <= path_bp < cum_bp[i+1].
        // Pure logic extracted to advance_step_cursor() for unit-testing
        // (see test_walker_*.) Records are sorted ascending in path_bp, so
        // this advances at most n_steps times in TOTAL across the slice.
        advance_step_cursor(path_bp, &cum_bp, &mut i, &mut cum);
        debug_assert!(cum_bp[i] <= path_bp);
        debug_assert!(path_bp < cum_bp[i + 1]);

        let curr_offset = path_bp - cum;

        // Dedup check (--dedup-read-node). After cursor advance, step_node_ids[i]
        // is the starting graph node and curr_offset is the intra-node offset
        // for this record's MEM walk. Two records sharing (read_id, read_st,
        // starting_node_id, offset) correspond to the same physical graph
        // position for the same read coordinate -- the intended duplicate
        // class emitted by find_mems --lightweight-tags (one record per tag
        // run sharing the same graph_pos within a MEM). Including offset
        // tightens dedup to true graph-position collisions; earlier versions
        // keyed on starting-node alone and silently collapsed records that
        // started at the same node but at different offsets.
        //
        // Key is packed into a u64 (vs. a 16-byte tuple); see walk_gfa() for
        // the bit layout and field-width asserts. HashSet<u64> hashes a
        // single word (cheaper than tuple hashing) and halves per-entry
        // memory: at HPRCv2 chr6 scale (60.7M unique triples), 625 MB vs.
        // 1.18 GB for the prior (u32, u32, usize) tuple key.
        if let Some(seen) = dedup_seen.as_mut() {
            // Pack (node_id - min_node_id) instead of raw node_id: see the
            // min_node_id parameter doc and walk_gfa()'s key-layout comment.
            // Subtraction is safe because step_node_ids[i] >= min_node_id by
            // construction (every step's node_id came from parse_gfa, which
            // also sourced min_id).
            let rel_node_id = step_node_ids[i] - min_node_id;
            let key = (r.read_id as u64)
                | ((r.read_st as u64) << 20)
                | ((rel_node_id as u64) << 29)
                | ((curr_offset as u64) << 54);
            // Layout uses all 64 bits; see walk_gfa() for field-width caps.
            if !seen.insert(key) {
                *dedup_skipped += 1;
                continue;
            }
        }

        path_str.clear();
        let path_len = traverse_nodes(&steps[i..], &step_node_ids[i..], callback, get_node_len, curr_offset, match_len, is_reverse, &mut path_str);
        let path_start = curr_offset;
        let path_end = curr_offset + match_len;
        if path_len < match_len {
            eprintln!("ERROR Path len {} << match_len {}. Path name: {}, Path str: {}", path_len, match_len, name, path_str);
            continue;
        }
        if path_len < path_end {
            eprintln!("ERROR: Path len {} <= path_end {}. Path name: {}, Path str: {}", path_len, path_end, name, path_str);
            continue;
        }

        total_gaf_entries += 1;

        if let Some(ref mut file) = gaf_output {
            write!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                   read_id, read_start, path_str, path_len, path_start, path_end, match_len, name)?;
        }
    }

    if path_bp_out_of_range > 0 {
        eprintln!("path {}: {} path_bp-out-of-range (of {} records)",
                  name, path_bp_out_of_range, num_matches);
    }
    if verbose {
        eprintln!("Processed all {} matches for path: {}\n", num_matches, name);
    }
    Ok((total_gaf_entries, 1))
}

fn walk_gfa(
    gfa_path: &str,
    path_pos_file: &str,
    path_to_seq_id_map: HashMap<String, usize>,
    path_starts: Vec<usize>,
    mut gaf_output: Option<BufWriter<File>>,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize,
    // When true, dedup (read_id, read_st, starting_node_id) triples across
    // the whole run. See process_path_matches() for details.
    dedup_read_node: bool,
    // Minimum raw GFA node_id (from parse_gfa). Forwarded to
    // process_path_matches so the dedup key can store a min_id-relative
    // node id and survive sparse-ID graphs (e.g. cactus-mc).
    min_node_id: usize,
    verbose: bool) -> std::io::Result<()>
{
    if let Some(ref mut file) = gaf_output {
        let header_line = "read_id\tread_st\tpath_str\tpath_len\tpath_st\tpath_end\tmatch_len\tpath_name";
        writeln!(file, "{}", header_line)?;
    }
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;

    // v2: _path_pos_v2.bin is a contiguous array of 16-byte Records. Cast in
    // place -- bytemuck guarantees alignment + layout via Pod/Zeroable.
    let records_load_start = Instant::now();
    let bytes = std::fs::read(path_pos_file)?;
    if bytes.len() % std::mem::size_of::<Record>() != 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("path_pos file size {} is not a multiple of {} (Record size). \
                     Expected v2 format (_path_pos_v2.bin). Is this a v1 file?",
                    bytes.len(), std::mem::size_of::<Record>())));
    }
    let records: &[Record] = bytemuck::cast_slice(&bytes);
    eprintln!("Loaded {} path_pos records ({} bytes, {} B/record)",
              records.len(), bytes.len(), std::mem::size_of::<Record>());
    eprintln!("PHASE_TIMING records_load_s: {:.6}", records_load_start.elapsed().as_secs_f64());

    let path_walk_start = Instant::now();
    let mut line = String::new();
    let mut total_gaf_entries = 0;
    let mut total_passes: u64 = 0;
    let mut max_passes: u64 = 0;
    let mut total_step_visits: u64 = 0;
    let mut nonempty_seq_ids: u64 = 0;

    // Dedup state. None disables dedup (default). Some(set) is shared across
    // ALL paths/seq_ids so cross-bucket duplicates (the common case from
    // lightweight find_mems) are caught. Sized to records.len() / 2 as an
    // initial guess: after dedup, count is typically 60-80% of input.
    //
    // Packed u64 key layout (low → high; all 64 bits used):
    //   bits [ 0..20):  read_id            (max 1,048,575 -- pipeline runs 500K reads)
    //   bits [20..29):  read_st            (max     511   -- short reads, ≤300 bp)
    //   bits [29..54):  node_id - min_id   (max 33,554,431 -- relative; see below)
    //   bits [54..64):  offset             (max     1023  -- node lengths capped ≤1024)
    //
    // The node-id field is stored RELATIVE to min_node_id (the smallest
    // segment id observed in parse_gfa). Raw GFA ids in some graph builds
    // (e.g. cactus-mc HPRCv2 chr6: min_id=162,720,608, max_id=171,559,458)
    // overflow 25 bits even though the actual id range fits in 24. Storing
    // (id - min_id) keeps the 25-bit field correct for both dense graphs
    // (min_id ≈ 1, pggb) and sparse-base graphs (min_id ≫ 1, cactus-mc).
    // The relative-id is identity-preserving on the dedup key: two records
    // hashing equal still represent the same physical (read_id, read_st,
    // node, offset) tuple.
    //
    // Why packed: HashSet<u64> is 8 B/slot vs. 16 B for (u32,u32,usize) and
    // hashes a single word. At 60.7M unique keys (HPRCv2 chr6), that's
    // ~625 MB vs. ~1.18 GB -- a 12% RSS shave on the 4.6 GB peak.
    //
    // Field-width asserts (below) convert silent overflow into loud failure.
    // Strand is intentionally NOT in the key: a fwd-bucket and rev-bucket
    // record landing on the same (node, offset) for the same (read, read_st)
    // would represent the same physical graph hit; if they did differ in
    // strand interpretation, find_mems should not have emitted both for a
    // single read coordinate.
    if dedup_read_node {
        let max_read_id = records.iter().map(|r| r.read_id).max().unwrap_or(0);
        let max_read_st = records.iter().map(|r| r.read_st).max().unwrap_or(0);
        assert!(max_read_id < (1u32 << 20),
            "read_id {} exceeds 20-bit packed dedup field (cap 1,048,575); \
             pipeline contract is 500K reads/run", max_read_id);
        assert!(max_read_st < (1u32 << 9),
            "read_st {} exceeds 9-bit packed dedup field (cap 511); \
             contract is short reads ≤300 bp", max_read_st);
    }
    let mut dedup_seen: Option<HashSet<u64>> = if dedup_read_node {
        eprintln!("Dedup mode: (read_id, read_st, starting_node_id, offset) — enabled");
        Some(HashSet::with_capacity(records.len() / 2 + 16))
    } else {
        None
    };
    let mut dedup_skipped: u64 = 0;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();

        // only parse paths
        if !line_str.starts_with('P') {
            continue;
        }

        // Parse segment line format: P<tab>p_name<tab>steps
        let mut fields = line_str.split('\t');

        let Some((name, steps)) = fields.next().and_then(|_type| {
            let name = fields.next()?;
            let steps = fields.next()?;
            Some((name, steps))
        }) else {
            eprintln!("Unable to parse GFA path: {}\n", line_str);
            continue;
        };

        // check if path name in map
        if !path_to_seq_id_map.contains_key(name) {
            eprintln!("Skipping path: {} not found\n", name);
            continue;
        }

        let positive_strand_seq_id = path_to_seq_id_map[name] * 2;  // Paths are represented by positive strands
        let steps: Vec<&str> = steps.split(',').collect();

        match process_path_matches(
            &steps,
            positive_strand_seq_id,
            name,
            false,
            &path_starts,
            records,
            &mut gaf_output,
            &mut callback,
            &mut get_node_len,
            &mut dedup_seen,
            &mut dedup_skipped,
            min_node_id,
            verbose,
        ) {
            Ok((entries, passes)) => {
                total_gaf_entries += entries;
                if passes > 0 {
                    nonempty_seq_ids += 1;
                    total_passes += passes;
                    max_passes = max_passes.max(passes);
                    total_step_visits += passes * steps.len() as u64;
                    if verbose {
                        eprintln!("PASSES\t{}\t{}\t{}\t{}", positive_strand_seq_id, steps.len(), entries, passes);
                    }
                }
            }
            Err(e) => eprintln!("Error processing path {}: {}", name, e),
        }

        let reverse_strand_seq_id = positive_strand_seq_id + 1;  // Paths are represented by positive strands

        let reverse_steps: Vec<&str> = steps.iter().rev().copied().collect();
        match process_path_matches(
            &reverse_steps,
            reverse_strand_seq_id,
            &format!("{}_reverse", name),
            true,
            &path_starts,
            records,
            &mut gaf_output,
            &mut callback,
            &mut get_node_len,
            &mut dedup_seen,
            &mut dedup_skipped,
            min_node_id,
            verbose,
        ) {
            Ok((entries, passes)) => {
                total_gaf_entries += entries;
                if passes > 0 {
                    nonempty_seq_ids += 1;
                    total_passes += passes;
                    max_passes = max_passes.max(passes);
                    total_step_visits += passes * steps.len() as u64;
                    if verbose {
                        eprintln!("PASSES\t{}\t{}\t{}\t{}", reverse_strand_seq_id, steps.len(), entries, passes);
                    }
                }
            }
            Err(e) => eprintln!("Error processing path {}: {}", name, e),
        }

    }
    if let Some(ref mut w) = gaf_output {
        w.flush()?;
    }
    eprintln!("PHASE_TIMING path_walk_s: {:.6}", path_walk_start.elapsed().as_secs_f64());
    eprintln!("---------------------");
    eprintln!("Total GAF entries: {}", total_gaf_entries);
    eprintln!("Path-scan passes: total={} over {} seq_ids (mean={:.1}, max={}); step-visits≈{}",
              total_passes, nonempty_seq_ids,
              total_passes as f64 / nonempty_seq_ids.max(1) as f64,
              max_passes, total_step_visits);
    if let Some(seen) = dedup_seen.as_ref() {
        let total_records = records.len();
        let kept = total_records as i64 - dedup_skipped as i64;
        eprintln!("Dedup: {} unique (read_id, read_st, starting_node_id, offset) keys; \
                   {} duplicates skipped ({:.2}% of {} records)",
                  seen.len(), dedup_skipped,
                  100.0 * dedup_skipped as f64 / total_records as f64,
                  total_records);
        // Sanity: number of inserts == number of kept records.
        // (kept == seen.len() if every kept record produced a new triple,
        //  which it does by construction.)
        if seen.len() as i64 != kept {
            eprintln!("WARN: dedup set size ({}) != kept records ({}); \
                       this should not happen — please report",
                      seen.len(), kept);
        }
    }
    Ok(())
}

/// Reads a file containing path names and creates a map from path name to line index
///
/// # Arguments
/// * `filename` - Path to the file containing path names (one per line)
///
/// # File Format
/// Each line contains a single path name.
/// Example:
/// ```
/// path1
/// path2
/// path3
/// ```
/// Results in: {"path1": 0, "path2": 1, "path3": 2}
///
/// # Returns
/// A HashMap where the key is the path name and value is the 0-based line index
fn read_path_name_indices(filename: &str, verbose: bool) -> std::io::Result<HashMap<String, usize>> {
    let mut path_indices = HashMap::new();
    let mut line_index = 0;

    for_each_line_in_file(filename, |line: &str| {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            return;
        }
        path_indices.insert(line.to_string(), line_index);
        line_index += 1;
    });

    if verbose {
        eprintln!("Path indices map:");
        for (path_name, index) in &path_indices {
            eprintln!("  {}: {}", path_name, index);
        }
    }
    eprintln!("Total paths: {}", path_indices.len());

    Ok(path_indices)
}


/// Project a GAF alignment file into coverage over GFA graph nodes
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Input GFA pangenome graph file (supports .gz/.bgz compression)
    #[arg(long)]
    gfa: String,
    /// Input GAF alignment file
    #[arg(short, long)]
    gaf: Option<String>,
    /// Scale coverage values by node length
    #[arg(short, long)]
    len_scale: bool,
    /// Emit graph coverage vector in a single column
    #[arg(short, long)]
    coverage_column: bool,
    /// Weight coverage by query group occurrences
    #[arg(short = 'w', long)]
    weight_queries: bool,
    /// Path positions binary file (find_mems v2 output: _path_pos_v2.bin).
    /// Records are 16 bytes each: path_bp, match_len_rev, read_st, read_id.
    #[arg(long)]
    path_pos: Option<String>,
    /// Sequence ID starts file
    #[arg(long)]
    seq_id_starts: Option<String>,
    /// Path names
    #[arg(long)]
    path_names: Option<String>,
    /// File prefix for output GAF file
    #[arg(long)]
    gaf_file_prefix: Option<String>,
    /// File prefix for coverage output (used when not generating GAF)
    #[arg(long)]
    coverage_prefix: Option<String>,
    /// Print per-path / per-record progress to stderr
    #[arg(short, long)]
    verbose: bool,
    /// Dedup (read_id, read_st, starting_node_id, offset) tuples when
    /// processing _path_pos_v2.bin input. Each unique 4-tuple contributes
    /// to coverage and GAF output exactly once; later records with the same
    /// key are skipped. Intended for use with find_mems --lightweight-tags,
    /// whose output contains intra-MEM same-graph-position duplicates (one
    /// per tag run sharing the same graph_pos within a MEM). Off by default:
    /// changes coverage semantics for any caller relying on per-haplotype
    /// duplicate records counting multiple times.
    ///
    /// Key is packed into a u64 with field-width caps enforced at load time:
    /// read_id < 2^20 (1M reads), read_st < 2^9 (512 bp), node_id < 2^25
    /// (33M nodes), offset < 2^10 (node length ≤1024 bp). Violations abort.
    #[arg(long)]
    dedup_read_node: bool,
    /// Worker threads for the path-walk phase (coverage-only mode).
    /// 1 = the original single-threaded walker. Output is byte-identical at
    /// any thread count: workers only compute coverage contributions, and a
    /// single-threaded merge applies them in GFA/bucket order, so dedup
    /// winners never depend on scheduling. --gaf falls back to the
    /// sequential walker regardless of this setting.
    #[arg(long, default_value_t = 8)]
    threads: usize,
}

// ============================================================================
// Parallel walker (coverage + dedup path)
// ============================================================================
//
// Why: PHASE_TIMING shows path_walk_s is ~93% of gafpack wall on HPRCv2-scale
// graphs (chr1: 177 s of 190 s). The dominant term inside it is the per-bucket
// rebuild of step_node_ids/cum_bp in process_path_matches -- ~5.97e9
// step-visits on chr1, and essentially INDEPENDENT of record count (mc-chr6
// shows 4,530,711,133 vs 4,530,698,953 step-visits for 1.74M vs 5.03M
// records). That work is per-path and embarrassingly parallel.
//
// Design: workers do the expensive per-bucket scan and return the coverage
// CONTRIBUTIONS they would have made; a single-threaded merge applies them in
// bucket order. The merge costs one hash op plus a few adds per record (~1-2 s
// against a 177 s scan), and buys three things:
//   * determinism -- output does not depend on thread scheduling, so coverage
//     MD5 remains a usable correctness gate;
//   * byte-identity with the sequential walker -- the merge visits buckets in
//     GFA order and records in index order, so dedup winners are the same
//     records the sequential walker would have kept. This matters because
//     ~1.05% of dedup keys on PGGB graphs have >=2 DISTINCT downstream walks
//     (measured), so a different winner yields different per-node coverage;
//   * no float-order question -- coverage accumulation stays single-threaded.
//
// Semantics preserved from process_path_matches (do not "simplify" these):
//   * the dedup insert happens BEFORE traversal, so a record that later fails
//     the path_len sanity checks still consumes its key;
//   * coverage IS applied for records that fail those checks -- traverse runs
//     first and the `continue` only skips the GAF-entry count;
//   * the path_len ERROR messages are emitted only for records that actually
//     reach traversal, i.e. dedup winners. Workers therefore record the error
//     condition and the merge decides whether to print it.
//
// --gaf is deliberately NOT supported here and falls back to the sequential
// walker: GAF output is a one-time validation artifact, while query.sh's
// production path is --coverage-prefix.

/// Per-record result from a worker. Coverage contributions live in the
/// bucket's flat `contribs` arena; this holds the range plus what the merge
/// needs to reproduce sequential stderr/counters.
#[derive(Clone, Copy)]
struct RecOut {
    key: u64,
    contrib_start: u32,
    contrib_len: u32,
    /// 0 = ok, 1 = path_len < match_len, 2 = path_len < path_end
    err: u8,
    path_len: u32,
    path_end: u32,
    match_len: u32,
}

/// One (path, orientation) bucket's worth of work, computed off-thread.
struct BucketOut {
    name: String,
    seq_id: usize,
    recs: Vec<RecOut>,
    /// (node_id - min_id, coverage_bp) pairs, indexed by RecOut ranges.
    contribs: Vec<(u32, u32)>,
    oob: u64,
    oob_msgs: Vec<String>,
    num_matches: usize,
    n_steps: usize,
}

/// traverse_nodes without the GAF path_str construction: pushes (node, bp)
/// contributions instead of invoking a callback. Iteration bound matches the
/// original's steps.zip(node_ids) -- both slices are the same length here.
fn traverse_collect(
    node_ids: &[usize],
    segment_lengths: &[usize],
    min_id: usize,
    start_offset: usize,
    match_len: usize,
    out: &mut Vec<(u32, u32)>,
) -> usize {
    if node_ids.is_empty() || match_len == 0 {
        return 0;
    }
    let mut total_path_length = 0usize;
    let mut remaining_len = match_len;
    for (i, &node_id) in node_ids.iter().enumerate() {
        let node_length = segment_lengths[node_id - min_id];
        total_path_length += node_length;
        let coverage = if i == 0 {
            min(remaining_len, node_length - start_offset)
        } else {
            min(remaining_len, node_length)
        };
        remaining_len -= coverage;
        out.push(((node_id - min_id) as u32, coverage as u32));
        if remaining_len == 0 {
            break;
        }
    }
    total_path_length
}


/// Build (node_ids, cum_bp) for a path in FORWARD order by parsing its steps.
/// This is the hot loop of the whole program: one integer parse plus one random
/// probe into segment_lengths per step, and step-visits reaches 5.97e9 on
/// HPRCv2.1 MC chr1. It must run once per path, never once per bucket.
#[inline]
fn build_forward_arrays(
    steps: &[&str],
    segment_lengths: &[usize],
    min_id: usize,
) -> (Vec<usize>, Vec<usize>) {
    let n = steps.len();
    let mut node_ids: Vec<usize> = Vec::with_capacity(n);
    let mut cum: Vec<usize> = Vec::with_capacity(n + 1);
    cum.push(0);
    for step in steps.iter() {
        let (seg, _orient) = step.split_at(step.len() - 1);
        let seg_id = seg.parse::<usize>().unwrap();
        node_ids.push(seg_id);
        cum.push(cum.last().unwrap() + segment_lengths[seg_id - min_id]);
    }
    (node_ids, cum)
}

/// Derive the REVERSE bucket's arrays from the forward ones, with no parsing
/// and no segment_lengths probes.
///
/// The reverse traversal visits nodes n-1, n-2, ... 0, so:
///     node_ids_rev[i] = node_ids_fwd[n-1-i]
///     cum_bp_rev[i]   = total - cum_bp_fwd[n-i]
/// Endpoints check out: cum_rev[0] = total - cum_fwd[n] = 0, and
/// cum_rev[n] = total - cum_fwd[0] = total. A node's length does not depend on
/// which direction it is traversed, which is what makes this exact rather than
/// approximate. Verified against a literal re-parse of the reversed step slice
/// in tests::reverse_arrays_match_reparse.
#[inline]
fn derive_reverse_arrays(
    node_ids_fwd: &[usize],
    cum_fwd: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let n = node_ids_fwd.len();
    let total = cum_fwd[n];
    let node_ids_rev: Vec<usize> = node_ids_fwd.iter().rev().copied().collect();
    let cum_rev: Vec<usize> = (0..=n).map(|i| total - cum_fwd[n - i]).collect();
    (node_ids_rev, cum_rev)
}

/// Pure per-bucket worker. No shared state, no I/O, no stderr.
#[allow(clippy::too_many_arguments)]
fn process_bucket(
    step_node_ids: &[usize],
    cum_bp: &[usize],
    seq_id: usize,
    name: &str,
    path_starts: &[usize],
    records: &[Record],
    segment_lengths: &[usize],
    min_id: usize,
    dedup: bool,
) -> BucketOut {
    debug_assert_eq!(cum_bp.len(), step_node_ids.len() + 1,
                     "cum_bp must have one more entry than step_node_ids");
    let st = path_starts[seq_id];
    let end = path_starts[seq_id + 1];
    let num_matches = end - st;

    let mut out = BucketOut {
        name: name.to_string(),
        seq_id,
        recs: Vec::new(),
        contribs: Vec::new(),
        oob: 0,
        oob_msgs: Vec::new(),
        num_matches,
        n_steps: step_node_ids.len(),
    };
    if num_matches == 0 {
        return out;
    }

    // Arrays are built once per PATH by the caller and shared between this
    // path's forward and reverse buckets -- see buckets_for_pline().
    let path_total_bp = *cum_bp.last().unwrap();

    out.recs.reserve(num_matches);
    let mut i: usize = 0;
    let mut cum: usize = 0;

    for idx in st..end {
        let r = &records[idx];
        let path_bp = r.path_bp as usize;
        let match_len = r.match_len as usize;

        if path_bp >= path_total_bp {
            out.oob += 1;
            if out.oob <= 5 {
                out.oob_msgs.push(format!(
                    "ERROR: path_bp {} >= path length {} on path {} (seq_id {}); skipping record idx {}",
                    path_bp, path_total_bp, name, seq_id, idx));
            }
            continue;
        }

        advance_step_cursor(path_bp, cum_bp, &mut i, &mut cum);
        debug_assert!(cum_bp[i] <= path_bp);
        debug_assert!(path_bp < cum_bp[i + 1]);
        let curr_offset = path_bp - cum;

        // Same packed layout as the sequential walker; see walk_gfa().
        let key = if dedup {
            let rel_node_id = step_node_ids[i] - min_id;
            (r.read_id as u64)
                | ((r.read_st as u64) << 20)
                | ((rel_node_id as u64) << 29)
                | ((curr_offset as u64) << 54)
        } else {
            0
        };

        let contrib_start = out.contribs.len() as u32;
        let path_len = traverse_collect(
            &step_node_ids[i..], segment_lengths, min_id,
            curr_offset, match_len, &mut out.contribs);
        let contrib_len = out.contribs.len() as u32 - contrib_start;
        let path_end = curr_offset + match_len;
        let err = if path_len < match_len { 1 } else if path_len < path_end { 2 } else { 0 };

        out.recs.push(RecOut {
            key, contrib_start, contrib_len, err,
            path_len: path_len as u32, path_end: path_end as u32,
            match_len: match_len as u32,
        });
    }
    out
}

/// Parse one GFA P-line into its two (forward, reverse) buckets. Returns an
/// empty vec for unparseable or unmapped paths, matching walk_gfa's behaviour
/// of warning and continuing.
fn buckets_for_pline(
    line_b: &[u8],
    path_to_seq_id_map: &HashMap<String, usize>,
    path_starts: &[usize],
    records: &[Record],
    segment_lengths: &[usize],
    min_id: usize,
    dedup: bool,
) -> Vec<BucketOut> {
    // Validation happens here, on a worker thread, rather than in the serial
    // read loop: read_line would UTF-8-validate all 26.4 GB single-threaded.
    let line_str = match std::str::from_utf8(line_b) {
        Ok(v) => v.trim(),
        Err(e) => {
            eprintln!("Skipping non-UTF-8 GFA path line: {}\n", e);
            return Vec::new();
        }
    };
    let mut fields = line_str.split('\t');
    let Some((name, steps_str)) = fields.next().and_then(|_type| {
        let name = fields.next()?;
        let steps = fields.next()?;
        Some((name, steps))
    }) else {
        eprintln!("Unable to parse GFA path: {}\n", line_str);
        return Vec::new();
    };
    let Some(&path_idx) = path_to_seq_id_map.get(name) else {
        eprintln!("Skipping path: {} not found\n", name);
        return Vec::new();
    };

    let seq_id = path_idx * 2;
    let rev_name = format!("{}_reverse", name);

    // Both buckets empty: skip the parse entirely. process_bucket used to
    // return before building its arrays in this case, so parsing here
    // unconditionally would ADD work on record-free paths (1,541 of the 9,634
    // buckets on mc-chr1). n_steps is only read by the merge when
    // num_matches > 0, so leaving it 0 here is sound.
    if path_starts[seq_id + 1] == path_starts[seq_id]
        && path_starts[seq_id + 2] == path_starts[seq_id + 1]
    {
        return vec![
            empty_bucket(name, seq_id),
            empty_bucket(&rev_name, seq_id + 1),
        ];
    }

    let steps: Vec<&str> = steps_str.split(',').collect();

    // Parse ONCE per path. The reverse bucket's arrays are pure arithmetic on
    // these, so the expensive per-step work (integer parse + a random probe
    // into a segment_lengths array larger than L3) is not repeated. `steps`
    // is dead after this: traverse_collect walks node_ids only, because the
    // parallel path never builds a GAF path_str.
    let (node_ids_fwd, cum_fwd) = build_forward_arrays(&steps, segment_lengths, min_id);
    drop(steps);

    let fwd = process_bucket(&node_ids_fwd, &cum_fwd, seq_id, name,
                             path_starts, records, segment_lengths, min_id, dedup);

    let (node_ids_rev, cum_rev) = derive_reverse_arrays(&node_ids_fwd, &cum_fwd);
    drop(node_ids_fwd);
    drop(cum_fwd);

    let rev = process_bucket(&node_ids_rev, &cum_rev, seq_id + 1, &rev_name,
                             path_starts, records, segment_lengths, min_id, dedup);
    vec![fwd, rev]
}

/// A bucket with no records: nothing for the merge to apply.
fn empty_bucket(name: &str, seq_id: usize) -> BucketOut {
    BucketOut {
        name: name.to_string(),
        seq_id,
        recs: Vec::new(),
        contribs: Vec::new(),
        oob: 0,
        oob_msgs: Vec::new(),
        num_matches: 0,
        n_steps: 0,
    }
}

/// Running state for the sequential merge.
struct MergeState {
    dedup_seen: Option<HashSet<u64>>,
    dedup_skipped: u64,
    total_gaf_entries: usize,
    total_passes: u64,
    max_passes: u64,
    total_step_visits: u64,
    nonempty_seq_ids: u64,
}

/// Apply one bucket's contributions. Visiting buckets in GFA order and
/// records in index order reproduces the sequential walker's dedup winners
/// exactly -- see the module comment above for why that matters.
fn merge_bucket(ms: &mut MergeState, coverage: &mut [f64], b: &BucketOut) {
    for m in &b.oob_msgs {
        eprintln!("{}", m);
    }
    if b.num_matches == 0 {
        return;
    }
    ms.nonempty_seq_ids += 1;
    ms.total_passes += 1;
    ms.max_passes = ms.max_passes.max(1);
    ms.total_step_visits += b.n_steps as u64;

    for rec in &b.recs {
        if let Some(seen) = ms.dedup_seen.as_mut() {
            if !seen.insert(rec.key) {
                ms.dedup_skipped += 1;
                continue;
            }
        }
        // Coverage is applied BEFORE the path_len checks, matching the
        // sequential walker (traverse_nodes runs first there too).
        let s = rec.contrib_start as usize;
        let e = s + rec.contrib_len as usize;
        for &(node, cov) in &b.contribs[s..e] {
            coverage[node as usize] += cov as f64;
        }
        match rec.err {
            1 => {
                eprintln!("ERROR Path len {} << match_len {}. Path name: {} (path_str omitted in parallel mode)",
                          rec.path_len, rec.match_len, b.name);
                continue;
            }
            2 => {
                eprintln!("ERROR: Path len {} <= path_end {}. Path name: {} (path_str omitted in parallel mode)",
                          rec.path_len, rec.path_end, b.name);
                continue;
            }
            _ => {}
        }
        ms.total_gaf_entries += 1;
    }

    if b.oob > 0 {
        eprintln!("path {}: {} path_bp-out-of-range (of {} records)",
                  b.name, b.oob, b.num_matches);
    }
}

#[allow(clippy::too_many_arguments)]
/// Fused single-pass parse + walk.
///
/// The sequential path reads the GFA twice: parse_gfa scans for S-lines, then
/// the walker scans again for P-lines. On HPRCv2.1 MC chr1 that is 2 x 26.4 GB,
/// and the parse leg alone measured ~9.9s of a 44.6s run.
///
/// GFA orders records H, S, L, P (verified on both chr1 graphs: S starts at
/// line 2, L at 11.8M, P at 28.1M), so one pass suffices: accumulate segment
/// lengths until the first P-line, finalise the dense tables there, then walk.
/// An S-line after finalisation would mean an incomplete node-length table, so
/// that aborts loudly rather than silently mis-covering.
///
/// Returns (coverage, segment_lengths, min_id) -- the caller needs the latter
/// two to write the coverage CSV.
fn parse_and_walk_parallel(
    gfa_path: &str,
    path_pos_file: &str,
    path_to_seq_id_map: &HashMap<String, usize>,
    path_starts: &[usize],
    dedup_read_node: bool,
    threads: usize,
    verbose: bool,
) -> std::io::Result<(Vec<f64>, Vec<usize>, usize)> {
    let records_load_start = Instant::now();
    let bytes = std::fs::read(path_pos_file)?;
    if bytes.len() % std::mem::size_of::<Record>() != 0 {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("path_pos file size {} is not a multiple of {} (Record size). \
                     Expected v2 format (_path_pos_v2.bin). Is this a v1 file?",
                    bytes.len(), std::mem::size_of::<Record>())));
    }
    let records: &[Record] = bytemuck::cast_slice(&bytes);
    eprintln!("Loaded {} path_pos records ({} bytes, {} B/record)",
              records.len(), bytes.len(), std::mem::size_of::<Record>());
    eprintln!("PHASE_TIMING records_load_s: {:.6}", records_load_start.elapsed().as_secs_f64());

    // Identical field-width contract to the sequential walker.
    if dedup_read_node {
        let max_read_id = records.iter().map(|r| r.read_id).max().unwrap_or(0);
        let max_read_st = records.iter().map(|r| r.read_st).max().unwrap_or(0);
        assert!(max_read_id < (1u32 << 20),
            "read_id {} exceeds 20-bit packed dedup field (cap 1,048,575); \
             pipeline contract is 500K reads/run", max_read_id);
        assert!(max_read_st < (1u32 << 9),
            "read_st {} exceeds 9-bit packed dedup field (cap 511); \
             contract is short reads <=300 bp", max_read_st);
    }

    // S-line accumulation state; finalised at the first P-line.
    let mut segments: Vec<(usize, usize)> = Vec::new();
    let mut seg_min_id = usize::MAX;
    let mut seg_max_id = 0usize;
    let mut segment_lengths: Vec<usize> = Vec::new();
    let mut coverage: Vec<f64> = Vec::new();
    let mut min_id: usize = 0;
    let mut finalized = false;
    let gfa_parse_start = Instant::now();

    let mut ms = MergeState {
        dedup_seen: if dedup_read_node {
            eprintln!("Dedup mode: (read_id, read_st, starting_node_id, offset) — enabled");
            Some(HashSet::with_capacity(records.len() / 2 + 16))
        } else { None },
        dedup_skipped: 0,
        total_gaf_entries: 0,
        total_passes: 0,
        max_passes: 0,
        total_step_visits: 0,
        nonempty_seq_ids: 0,
    };

    let path_walk_start = Instant::now();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("failed to build rayon thread pool");
    // Batch is deliberately small: peak memory is one batch of step arrays
    // (~12 MB per HPRC chr1 bucket) plus its contributions, not the whole GFA.
    let batch_size = (threads * 2).max(2);
    eprintln!("Parallel walker: {} threads, batch {} P-lines", threads, batch_size);

    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    // Raw bytes, not read_line-into-String: the walk phase rescans the whole
    // GFA (28.08M lines / 26.4 GB on chr1) and only 4,817 of those lines are
    // P-lines. read_line would UTF-8-validate every byte on the serial path.
    //
    // P-lines are read DIRECTLY into a per-batch arena rather than into a
    // scratch buffer that is then cloned into the batch. read_until appends,
    // so we can speculatively read every line at the arena's tail and simply
    // truncate back when it turns out not to be a P-line. That removes an
    // entire copy of the ~25 GB of P-line payload per run; the earlier
    // scratch+clone shape paid for it twice.
    let mut arena: Vec<u8> = Vec::new();
    let mut spans: Vec<(usize, usize)> = Vec::with_capacity(batch_size);
    // Serial-phase accounting: merge_secs is the ordered-merge cost (the price
    // of determinism); read_secs is the single-threaded GFA line scan.
    // Serial-phase accounting. Measured by subtraction rather than per-line
    // Instant::now() calls: the scan touches ~28M lines on chr1, so two clock
    // reads per line would be ~1s of pure instrumentation. read_secs therefore
    // covers ALL line reading -- including the S/L lines that are scanned and
    // discarded, which an earlier version of this counter silently omitted by
    // accumulating after the `continue`.
    let mut merge_secs = 0.0f64;
    let mut dispatch_secs = 0.0f64;
    let mut lines_scanned: u64 = 0;
    let mut p_lines: u64 = 0;

    loop {
        let line_start = arena.len();
        let bytes_read = reader.read_until(b'\n', &mut arena)?;
        let eof = bytes_read == 0;
        if !eof {
            lines_scanned += 1;
            // GFA record type is the first byte of the line; no leading
            // whitespace is permitted by the spec.
            let rec_type = arena[line_start];
            if rec_type == b'S' {
                assert!(!finalized,
                    "GFA S-line found after the first P-line: the node-length \
                     table was already finalised, so this segment would be \
                     missing from coverage. Re-run with --threads 1 for this graph.");
                let line_b = trim_eol(&arena[line_start..]);
                let mut fields = line_b.split(|&c| c == b'\t');
                if let Some((id_b, seq)) = fields.next().and_then(|_t| {
                    let id_b = fields.next()?;
                    let seq = fields.next()?;
                    Some((id_b, seq))
                }) {
                    let id = parse_usize_bytes(id_b).expect("non-numeric GFA segment id");
                    seg_min_id = seg_min_id.min(id);
                    seg_max_id = seg_max_id.max(id);
                    segments.push((id, seq.len()));
                }
                arena.truncate(line_start);
                continue;
            }
            if rec_type != b'P' {
                arena.truncate(line_start);
                continue;
            }
            // First P-line: every S-line has been seen, so build dense tables.
            if !finalized {
                min_id = if seg_min_id == usize::MAX { 0 } else { seg_min_id };
                let num_segments = if segments.is_empty() { 0 } else { seg_max_id - min_id + 1 };
                segment_lengths = vec![0; num_segments];
                for (id, len) in segments.drain(..) {
                    segment_lengths[id - min_id] = len;
                }
                if dedup_read_node {
                    let max_rel_node_id = num_segments.saturating_sub(1);
                    let max_node_len = segment_lengths.iter().copied().max().unwrap_or(0);
                    assert!(max_rel_node_id < (1usize << 25),
                        "node-id range (num_segments={}) exceeds 25-bit packed dedup \
                         field (cap 33,554,431). min_id={}", num_segments, min_id);
                    assert!(max_node_len <= 1024,
                        "node length {} exceeds 10-bit packed offset field \
                         (offsets 0..1023 -> node length cap 1024)", max_node_len);
                }
                coverage = vec![0.0; num_segments];
                eprintln!("PHASE_TIMING gfa_parse_s: {:.6}   (fused into the walk pass)",
                          gfa_parse_start.elapsed().as_secs_f64());
                finalized = true;
            }
            spans.push((line_start, arena.len()));
            p_lines += 1;
        }
        if spans.len() >= batch_size || (eof && !spans.is_empty()) {
            let dispatch_start = Instant::now();
            let arena_ref: &[u8] = &arena;
            let outs: Vec<BucketOut> = pool.install(|| {
                spans.par_iter()
                    .flat_map_iter(|&(a, b)| buckets_for_pline(
                        &arena_ref[a..b], path_to_seq_id_map, path_starts, records,
                        &segment_lengths, min_id, dedup_read_node))
                    .collect()
            });
            dispatch_secs += dispatch_start.elapsed().as_secs_f64();
            let merge_start = Instant::now();
            for b in &outs {
                if verbose {
                    eprintln!("PASSES\t{}\t{}\t{}\t{}", b.seq_id, b.n_steps, b.recs.len(), 1);
                }
                merge_bucket(&mut ms, &mut coverage, b);
            }
            merge_secs += merge_start.elapsed().as_secs_f64();
            spans.clear();
            arena.clear();
        }
        if eof {
            break;
        }
    }

    eprintln!("PHASE_TIMING path_walk_s: {:.6}", path_walk_start.elapsed().as_secs_f64());
    let walk_total = path_walk_start.elapsed().as_secs_f64();
    eprintln!("PHASE_TIMING   walk_dispatch_s: {:.6}   (parallel per-bucket work)", dispatch_secs);
    eprintln!("PHASE_TIMING   walk_merge_s: {:.6}   (ordered merge = cost of determinism)", merge_secs);
    eprintln!("PHASE_TIMING   walk_read_s: {:.6}   (single-threaded GFA line scan, {} lines, {} of them P)",
              (walk_total - dispatch_secs - merge_secs).max(0.0), lines_scanned, p_lines);
    eprintln!("---------------------");
    eprintln!("Total GAF entries: {}", ms.total_gaf_entries);
    eprintln!("Path-scan passes: total={} over {} seq_ids (mean={:.1}, max={}); step-visits≈{}",
              ms.total_passes, ms.nonempty_seq_ids,
              ms.total_passes as f64 / ms.nonempty_seq_ids.max(1) as f64,
              ms.max_passes, ms.total_step_visits);
    if let Some(seen) = ms.dedup_seen.as_ref() {
        let total_records = records.len();
        let kept = total_records as i64 - ms.dedup_skipped as i64;
        eprintln!("Dedup: {} unique (read_id, read_st, starting_node_id, offset) keys; \
                   {} duplicates skipped ({:.2}% of {} records)",
                  seen.len(), ms.dedup_skipped,
                  100.0 * ms.dedup_skipped as f64 / total_records as f64,
                  total_records);
        if seen.len() as i64 != kept {
            eprintln!("WARN: dedup set size ({}) != kept records ({}); \
                       this should not happen — please report",
                      seen.len(), kept);
        }
    }
    if !finalized {
        // GFA with no P-lines at all: still hand back a consistent table.
        min_id = if seg_min_id == usize::MAX { 0 } else { seg_min_id };
        let num_segments = if segments.is_empty() { 0 } else { seg_max_id - min_id + 1 };
        segment_lengths = vec![0; num_segments];
        for (id, len) in segments.drain(..) {
            segment_lengths[id - min_id] = len;
        }
        coverage = vec![0.0; num_segments];
    }
    Ok((coverage, segment_lengths, min_id))
}

fn main() {
    let args = Args::parse();
    let gfa_file = &args.gfa;

    // Check if we're in CSV processing mode
    if let (Some(path_pos_file), Some(seq_id_starts_file)) = (&args.path_pos, &args.seq_id_starts) {
        // let path_pos_data = read_path_pos_csv(path_pos_file).unwrap();
        eprintln!("Processing gfa file: {}", gfa_file);
        eprintln!("Args: {}, {}", path_pos_file, seq_id_starts_file);
        let path_names_file_name = &args.path_names.unwrap();
        eprintln!("Path names file: {}", path_names_file_name);
        // let seq_id_starts_map = read_seq_id_starts_map(seq_id_starts_file).unwrap();

        let seq_starts = read_cumulative_starts(seq_id_starts_file).unwrap();
        let path_to_seq_id_map = read_path_name_indices(path_names_file_name, args.verbose).unwrap();

        // Create output file if file_prefix is provided
        let gaf_output = if let Some(prefix) = &args.gaf_file_prefix {
            let gaf_filename = format!("{}.gaf", prefix);
            let f = File::create(gaf_filename).expect("Could not create GAF output file");
            Some(BufWriter::with_capacity(1 << 20, f))
        } else {
            None
        };

        // Parallel walker handles the production path (coverage-only), and
        // fuses the S-line parse into the same GFA pass. --gaf stays on the
        // two-pass sequential walker: it is a one-time validation artifact and
        // keeping GAF line order identical there is worth more than the
        // speedup. --threads 1 also routes to sequential, giving a trivial A/B.
        let use_parallel = args.threads > 1 && gaf_output.is_none();

        let (coverage, segment_lengths, min_id) = if use_parallel {
            match parse_and_walk_parallel(
                &gfa_file, path_pos_file, &path_to_seq_id_map, &seq_starts,
                args.dedup_read_node, args.threads, args.verbose)
            {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("Error processing GFA file with path positions: {}", e);
                    std::process::exit(1);
                }
            }
        } else {
            // Two-pass: parse_gfa scans for S-lines, walk_gfa rescans for P.
            let gfa_parse_start = Instant::now();
            let (segment_lengths, min_id) = parse_gfa(&gfa_file).unwrap();
            eprintln!("PHASE_TIMING gfa_parse_s: {:.6}", gfa_parse_start.elapsed().as_secs_f64());
            let num_segments = segment_lengths.len();

            // Field-width contract for the packed dedup key (see walk_gfa()).
            // Checked here (not in parse_gfa) because the GAF-input mode
            // doesn't use the packed key and shouldn't be constrained by it.
            //
            // The node-id field stores (node_id - min_id), so the constraint
            // is on the id RANGE (num_segments), not the absolute max id.
            // This lets graphs with sparse high-base ids (cactus-mc:
            // min_id≈163M) work as long as the id span fits 25 bits.
            if args.dedup_read_node {
                let max_rel_node_id = num_segments - 1;
                let max_node_len = segment_lengths.iter().copied().max().unwrap_or(0);
                assert!(max_rel_node_id < (1usize << 25),
                    "node-id range (num_segments={}) exceeds 25-bit packed \
                     dedup field (cap 33,554,431). min_id={}, max_id={}.",
                    num_segments, min_id, min_id + num_segments - 1);
                assert!(max_node_len <= 1024,
                    "node length {} exceeds 10-bit packed offset field \
                     (offsets 0..1023 → node length cap 1024)",
                    max_node_len);
            }

            let mut coverage: Vec<f64> = vec![0.0; num_segments];
            if let Err(e) = walk_gfa(&gfa_file, path_pos_file, path_to_seq_id_map, seq_starts, gaf_output, |node_id, len| {
                coverage[node_id - min_id] += len as f64;
            }, |node_id| segment_lengths[node_id - min_id], args.dedup_read_node, min_id, args.verbose) {
                eprintln!("Error processing GFA file with path positions: {}", e);
                std::process::exit(1);
            }
            (coverage, segment_lengths, min_id)
        };

        let coverage_write_start = Instant::now();
        let output_filename = format!("{}_coverage.csv",
                                      args.gaf_file_prefix.as_deref()
                                          .or(args.coverage_prefix.as_deref())
                                          .unwrap_or("output"));
        let mut output_file = BufWriter::with_capacity(
            1 << 20,
            File::create(&output_filename).expect("Could not create coverage output file"),
        );

        // Write header
        writeln!(output_file, "node_id,node_coverage").unwrap();

        // Write coverage data
        for (i, v) in coverage.into_iter().enumerate() {
            let node_id = min_id + i;
            let coverage_value = if args.len_scale {
                v / segment_lengths[i] as f64
            } else {
                v
            };
            writeln!(output_file, "{},{:.2}", node_id, coverage_value).unwrap();
        }
        output_file.flush().expect("flush coverage output");
        eprintln!("PHASE_TIMING coverage_write_s: {:.6}", coverage_write_start.elapsed().as_secs_f64());

        eprintln!("Coverage data written to: {}", output_filename);
        return;
        
        
        
        // return;
        // eprint!("#sample");
        // for n in min_id..min_id + num_segments {
        //     eprint!("\tnode.{}", n);
        // }
        // eprintln!();
        // for (i, v) in coverage.into_iter().enumerate() {
        //     print!(
        //         "\t{}",
        //         if args.len_scale {
        //             v / segment_lengths[i] as f64
        //         } else {
        //             v
        //         }
        //     );
        // }
        // eprintln!();
        // return;
    }

    // Original coverage calculation logic
    let gaf_file = args.gaf.expect("GAF file is required for coverage calculation");

    // Parse GFA file
    let (segment_lengths, min_id) = parse_gfa(&gfa_file).unwrap();
    let num_segments = segment_lengths.len();       // segment -> node

    let mut coverage: Vec<f64> = vec![0.0; num_segments];

    if args.weight_queries {
        // First pass: count query occurrences
        let mut query_counts: HashMap<String, usize> = HashMap::new();
        for_each_line_in_file(&gaf_file, |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            if fields.len() >= 4 {
                let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
                *query_counts.entry(query_key).or_insert(0) += 1;
            }
        });

        // Second pass: calculate coverage with query count adjustment
        for_each_line_in_file(&gaf_file, |l: &str| {
            let fields: Vec<&str> = l.split('\t').collect();
            let query_key = format!("{}:{}:{}", fields[0], fields[2], fields[3]);
            let count = query_counts.get(&query_key).unwrap_or(&1);

            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64 / *count as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        });
    } else {
        // Single pass without weighting
        for_each_line_in_file(&gaf_file, |l: &str| {
            for_each_step(
                l,
                |node_id, len| {
                    coverage[node_id - min_id] += len as f64;
                },
                |node_id| segment_lengths[node_id - min_id],
            );
        });
    }

    if args.coverage_column {
        println!("##sample: {}", gaf_file);
        println!("#coverage");
        for (i, v) in coverage.into_iter().enumerate() {
            eprintln!(
                "{}",
                if args.len_scale {
                    v / segment_lengths[i] as f64
                } else {
                    v
                }
            );
        }
    } else {
        print!("#sample");
        for n in min_id..min_id + num_segments {
            print!("\tnode.{}", n);
        }
        eprintln!();
        print!("{}", gaf_file);
        for (i, v) in coverage.into_iter().enumerate() {
            print!(
                "\t{}",
                if args.len_scale {
                    v / segment_lengths[i] as f64
                } else {
                    v
                }
            );
        }
        eprintln!();
    }
}

// =============================================================================
// v2 unit tests (see PLAN_find_mems_binary_io_v2.md, Gate 1).
// Run with:  cargo test
// =============================================================================
#[cfg(test)]
mod tests {
    use super::*;
    use bytemuck::cast_slice;

    // -------------------------------------------------------------------------
    // Gate 1a: record round-trip + bit-31 packing
    // -------------------------------------------------------------------------

    #[test]
    fn record_size_is_16_bytes() {
        assert_eq!(std::mem::size_of::<Record>(), 16);
    }

    #[test]
    fn record_decodes_known_bytes() {
        // Hand-encode 2 records, little-endian: path_bp, match_len, read_st, read_id
        // Record 0: path_bp=0x01020304, match_len=30,       read_st=7, read_id=42
        // Record 1: path_bp=u32::MAX,   match_len=20_000,   read_st=0, read_id=99
        let bytes: [u8; 32] = [
            // Record 0
            0x04, 0x03, 0x02, 0x01,             // path_bp = 0x01020304
            0x1E, 0x00, 0x00, 0x00,             // match_len = 30
            0x07, 0x00, 0x00, 0x00,             // read_st = 7
            0x2A, 0x00, 0x00, 0x00,             // read_id = 42
            // Record 1
            0xFF, 0xFF, 0xFF, 0xFF,             // path_bp = u32::MAX
            0x20, 0x4E, 0x00, 0x00,             // match_len = 20_000
            0x00, 0x00, 0x00, 0x00,             // read_st = 0
            0x63, 0x00, 0x00, 0x00,             // read_id = 99
        ];
        let records: &[Record] = cast_slice(&bytes);
        assert_eq!(records.len(), 2);

        assert_eq!(records[0].path_bp, 0x01020304);
        assert_eq!(records[0].match_len, 30);
        assert_eq!(records[0].read_st, 7);
        assert_eq!(records[0].read_id, 42);

        assert_eq!(records[1].path_bp, u32::MAX);
        assert_eq!(records[1].match_len, 20_000);
        assert_eq!(records[1].read_st, 0);
        assert_eq!(records[1].read_id, 99);
    }

    #[test]
    fn record_match_len_full_u32_range() {
        // Now that is_rev is gone, match_len uses the full u32 range.
        // Reads are ~20Kbp in practice but we shouldn't artificially cap.
        let cases: &[u32] = &[1, 20_000, 0x7FFF_FFFF, u32::MAX];
        for &mlen in cases {
            let r = Record { path_bp: 0, match_len: mlen, read_st: 0, read_id: 0 };
            assert_eq!(r.match_len, mlen, "match_len round-trip failed for {}", mlen);
        }
    }

    // -------------------------------------------------------------------------
    // Gate 1b: walker derivation edge cases (advance_step_cursor)
    //
    // Synthetic fixture:
    //   step_node_ids = [10, 20, 30, 40]
    //   node lengths  = [ 5,  7,  1,  3]
    //   cum_bp        = [ 0,  5, 12, 13, 16]   (path_total_bp = 16)
    //
    // For each test case we feed a single path_bp through the cursor and verify
    // the derived (step_index, local_offset). step_node_ids[i] is the node_id,
    // path_bp - cum_bp[i] is the offset. This is the exact derivation
    // process_path_matches does.
    // -------------------------------------------------------------------------


    // ---- reverse-array derivation -----------------------------------------
    // derive_reverse_arrays() replaces re-parsing a reversed step slice. The
    // reference below IS the old code path (reverse the slice, then parse), so
    // these assert exact equivalence rather than merely plausible behaviour.

    fn reference_reverse(steps: &[&str], lens: &[usize], min_id: usize)
        -> (Vec<usize>, Vec<usize>)
    {
        let rev: Vec<&str> = steps.iter().rev().copied().collect();
        build_forward_arrays(&rev, lens, min_id)
    }

    fn check_reverse(steps: &[&str], lens: &[usize], min_id: usize) {
        let (ids_f, cum_f) = build_forward_arrays(steps, lens, min_id);
        let (ids_d, cum_d) = derive_reverse_arrays(&ids_f, &cum_f);
        let (ids_r, cum_r) = reference_reverse(steps, lens, min_id);
        assert_eq!(ids_d, ids_r, "node_ids mismatch for {:?}", steps);
        assert_eq!(cum_d, cum_r, "cum_bp mismatch for {:?}", steps);
        // Structural invariants the walker relies on.
        assert_eq!(cum_d.len(), ids_d.len() + 1);
        assert_eq!(cum_d[0], 0);
        assert_eq!(*cum_d.last().unwrap(), *cum_f.last().unwrap());
        assert!(cum_d.windows(2).all(|w| w[0] <= w[1]), "cum_bp must be non-decreasing");
    }

    #[test]
    fn reverse_arrays_match_reparse() {
        // ids 5..=9, lengths indexed by (id - min_id)
        let lens = vec![3, 1, 7, 2, 5];
        check_reverse(&["5+", "7-", "9+", "6+", "8-"], &lens, 5);
    }

    #[test]
    fn reverse_arrays_single_step() {
        let lens = vec![11];
        check_reverse(&["0+"], &lens, 0);
        check_reverse(&["0-"], &lens, 0);
    }

    #[test]
    fn reverse_arrays_uniform_unit_lengths() {
        let lens = vec![1; 6];
        check_reverse(&["0+", "1+", "2+", "3+", "4+", "5+"], &lens, 0);
    }

    #[test]
    fn reverse_arrays_sparse_min_id() {
        // cactus-mc style: ids far from zero (observed min_id ~163M).
        let min_id = 162_720_608usize;
        let lens = vec![4, 9, 2, 6];
        let steps = ["162720609+", "162720611-", "162720608+", "162720610+"];
        check_reverse(&steps, &lens, min_id);
    }

    #[test]
    fn reverse_arrays_repeated_nodes() {
        // Same node visited several times, both orientations.
        let lens = vec![5, 8];
        check_reverse(&["0+", "1-", "0-", "0+", "1+"], &lens, 0);
    }

    #[test]
    fn reverse_arrays_randomized() {
        // Deterministic LCG; no dev-dependency needed.
        let mut state: u64 = 0x2026_08_31;
        let mut next = move || { state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (state >> 33) as usize };
        for _trial in 0..500 {
            let n_nodes = 1 + next() % 12;
            let lens: Vec<usize> = (0..n_nodes).map(|_| 1 + next() % 1024).collect();
            let n_steps = 1 + next() % 40;
            let owned: Vec<String> = (0..n_steps)
                .map(|_| format!("{}{}", next() % n_nodes, if next() % 2 == 0 { "+" } else { "-" }))
                .collect();
            let steps: Vec<&str> = owned.iter().map(|s| s.as_str()).collect();
            check_reverse(&steps, &lens, 0);
        }
    }

    fn fixture_cum_bp() -> Vec<usize> {
        vec![0, 5, 12, 13, 16]
    }
    fn fixture_node_ids() -> Vec<usize> {
        vec![10, 20, 30, 40]
    }
    const FIXTURE_PATH_TOTAL_BP: usize = 16;

    /// Run a single record through the cursor starting at i=0, return (i, offset)
    /// or None if rejected by the path_total_bp guard.
    fn derive_step(path_bp: usize, cum_bp: &[usize], path_total: usize) -> Option<(usize, usize)> {
        if path_bp >= path_total { return None; }
        let mut i = 0usize;
        let mut cum = 0usize;
        advance_step_cursor(path_bp, cum_bp, &mut i, &mut cum);
        Some((i, path_bp - cum))
    }

    #[test]
    fn walker_e1_first_bp_of_first_node() {
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(0, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (10, 0));
    }

    #[test]
    fn walker_e2_last_bp_of_first_node() {
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(4, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (10, 4));
    }

    #[test]
    fn walker_e3_exact_step_boundary_into_new_step() {
        // Critical off-by-one: path_bp == cum_bp[i+1] must land on step i+1,
        // not step i. The while-condition `path_bp >= cum_bp[i+1]` enforces this.
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(5, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (20, 0));
    }

    #[test]
    fn walker_e4_first_bp_inside_multi_bp_node() {
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(6, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (20, 1));
    }

    #[test]
    fn walker_e5_last_bp_of_multi_bp_node() {
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(11, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (20, 6));
    }

    #[test]
    fn walker_e6_entry_to_single_bp_node() {
        // Cursor must advance to step 2 (length-1 node 30).
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(12, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (30, 0));
    }

    #[test]
    fn walker_e7_exit_single_bp_node_into_next() {
        // Cursor must advance to step 3, crossing the 1-bp node 30 cleanly.
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(13, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (40, 0));
    }

    #[test]
    fn walker_e8_last_valid_bp_on_path() {
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let (i, off) = derive_step(15, &cum_bp, FIXTURE_PATH_TOTAL_BP).unwrap();
        assert_eq!((nodes[i], off), (40, 2));
    }

    #[test]
    fn walker_e9_one_past_end_rejected() {
        let cum_bp = fixture_cum_bp();
        assert!(derive_step(16, &cum_bp, FIXTURE_PATH_TOTAL_BP).is_none());
    }

    #[test]
    fn walker_e10_past_end_rejected() {
        let cum_bp = fixture_cum_bp();
        assert!(derive_step(17, &cum_bp, FIXTURE_PATH_TOTAL_BP).is_none());
    }

    #[test]
    fn walker_e11_u32_max_rejected_no_overflow() {
        // u32::MAX as usize must be safely rejected without arithmetic overflow.
        let cum_bp = fixture_cum_bp();
        assert!(derive_step(u32::MAX as usize, &cum_bp, FIXTURE_PATH_TOTAL_BP).is_none());
    }

    // -------- Monotonic-cursor stress: simulate a full slice of records --------

    #[test]
    fn walker_monotonic_cursor_across_slice() {
        // Records sorted by path_bp ascending, including repeats, exact
        // boundaries, gaps. Mirrors the real per-bucket loop in
        // process_path_matches. Asserts cursor is non-decreasing throughout
        // AND derives the correct (node_id, offset) per record.
        let cum_bp = fixture_cum_bp();
        let nodes = fixture_node_ids();
        let path_bps: &[usize] = &[0, 1, 1, 2, 5, 5, 5, 6, 12, 12, 13, 15];
        let expected: &[(usize, usize)] = &[
            (10, 0), (10, 1), (10, 1), (10, 2),
            (20, 0), (20, 0), (20, 0), (20, 1),
            (30, 0), (30, 0),
            (40, 0), (40, 2),
        ];
        assert_eq!(path_bps.len(), expected.len());

        let mut i = 0usize;
        let mut cum = 0usize;
        let mut last_i = 0usize;

        for (k, &pbp) in path_bps.iter().enumerate() {
            advance_step_cursor(pbp, &cum_bp, &mut i, &mut cum);
            assert!(i >= last_i,
                "cursor moved backwards at record {}: i={} last_i={} path_bp={}",
                k, i, last_i, pbp);
            assert_eq!(cum, cum_bp[i],
                "cum out of sync at record {}: cum={} cum_bp[{}]={}", k, cum, i, cum_bp[i]);
            let off = pbp - cum;
            assert_eq!((nodes[i], off), expected[k],
                "wrong derivation at record {} (path_bp={})", k, pbp);
            last_i = i;
        }
    }

    #[test]
    fn walker_isolation_between_calls() {
        // The forward and reverse passes share `records` but call
        // process_path_matches with DIFFERENT seq_ids => disjoint slices,
        // and each call gets its own freshly-zeroed cursor. Verify by
        // running advance_step_cursor twice "from scratch" on the same path,
        // each time starting from i=0/cum=0.
        let cum_bp = fixture_cum_bp();

        // Call 1: walk records [0, 5, 12, 15]
        let mut i = 0usize;
        let mut cum = 0usize;
        for &pbp in &[0usize, 5, 12, 15] {
            advance_step_cursor(pbp, &cum_bp, &mut i, &mut cum);
        }
        assert_eq!(i, 3);
        assert_eq!(cum, 13);

        // Call 2: simulate the second orientation -- cursor MUST start fresh
        let mut i = 0usize;
        let mut cum = 0usize;
        advance_step_cursor(0, &cum_bp, &mut i, &mut cum);
        assert_eq!((i, cum), (0, 0),
            "second call leaked state from first; cursor should restart at i=0");
    }

    #[test]
    fn walker_single_step_path() {
        // Degenerate path with 1 step, length 5. cum_bp = [0, 5].
        // path_bps 0..=4 must all map to step 0, offset 0..=4. bp 5 must reject.
        let cum_bp = vec![0usize, 5];
        for pbp in 0..5 {
            let mut i = 0usize;
            let mut cum = 0usize;
            advance_step_cursor(pbp, &cum_bp, &mut i, &mut cum);
            assert_eq!(i, 0, "single-step path: i should never advance");
            assert_eq!(pbp - cum, pbp);
        }
        assert!(derive_step(5, &cum_bp, 5).is_none());
    }
}
