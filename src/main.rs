use std::cmp::min;
use clap::Parser;
use flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

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

    if path
        .extension()
        .is_some_and(|ext| ext == "gz" || ext == "bgz")
    {
        let decoder = GzDecoder::new(file);
        let buf_reader = BufReader::new(decoder);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(file);
        Ok(Box::new(buf_reader))
    }
}

/// Parse GFA file and extract segment information
/// Returns (segment_lengths, min_id) where segment_lengths[id - min_id] gives the length
fn parse_gfa(gfa_path: &str) -> std::io::Result<(Vec<usize>, usize)> {
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;
    let mut line = String::new();
    let mut segments_map = HashMap::new();
    let mut min_id = usize::MAX;
    let mut max_id = 0;

    loop {
        line.clear();
        let bytes_read = reader.read_line(&mut line)?;
        if bytes_read == 0 {
            break;
        }

        let line_str = line.trim();

        // Only process segment lines
        if !line_str.starts_with('S') {
            continue;
        }

        // Parse segment line format: S<tab>id<tab>sequence
        let mut fields = line_str.split('\t');
        let Some((id_str, seq)) = fields.next().and_then(|_type| {
            let id_str = fields.next()?;
            let seq = fields.next()?;
            Some((id_str, seq))
        }) else {
            continue;
        };

        // Parse segment ID
        let id = id_str.parse::<usize>().unwrap();
        min_id = min_id.min(id);
        max_id = max_id.max(id);
        segments_map.insert(id, seq.len());
    }

    // Create a dense vector for O(1) access
    let num_segments = max_id - min_id + 1;
    let mut segment_lengths = vec![0; num_segments];

    for (id, len) in segments_map {
        segment_lengths[id - min_id] = len;
    }

    Ok((segment_lengths, min_id))
}


fn traverse_nodes(
    steps: &[&str],
    callback: &mut impl FnMut(usize, usize),
    get_node_len: &mut impl FnMut(usize) -> usize,
    start_offset: usize,
    alignment_len: usize,
) -> (String, usize)
{
    if steps.len() == 0 || alignment_len == 0 {
        return (String::new(), 0);
    }

    let mut total_path_length = 0;
    let mut traversed_nodes = Vec::new();
    let mut remaining_len = alignment_len;

    for (i, step) in steps.iter().enumerate() {
        let (seg, _orient) = step.split_at(step.len() - 1);
        let node_id = seg.parse::<usize>().unwrap();
        let node_length = get_node_len(node_id);
        total_path_length += node_length;
        traversed_nodes.push(node_id.to_string());

        let coverage = if i == 0 {
            // First node: account for start_offset
            min(remaining_len, node_length - start_offset)
        } else {
            // Subsequent nodes: use full node length or remaining length
            min(remaining_len, node_length)
        };

        remaining_len -= coverage;
        callback(node_id, node_length);

        if remaining_len <= 0 {
            break;
        }
    }

    (traversed_nodes.join(","), total_path_length)
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

        // println!("Line: {}", line);

        if line.is_empty() {
            continue;
        }
        // println!("Line bytes: {:?}", line.as_bytes());

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

        println!("Storing: seq_id: {}, val: {}, {}", seq_id, cumulative_start, count);

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

    eprintln!("Cumulative starts vector:");
    for (index, value) in cumulative_starts.iter().enumerate() {
        eprintln!("  [{}]: {}", index, value);
    }
    eprintln!("Total entries: {}", cumulative_starts.len());

    Ok(cumulative_starts)
}

fn walk_gfa(
    gfa_path: &str,
    path_pos_file: &str,
    path_to_seq_id_map: HashMap<String, usize>,
    path_starts: Vec<usize>,
    mut callback: impl FnMut(usize, usize),
    mut get_node_len: impl FnMut(usize) -> usize) -> std::io::Result<()>
{
    let path = Path::new(gfa_path);
    let mut reader = create_reader(path)?;

    let path_pos_path = Path::new(path_pos_file);

    let mut line = String::new();

    // dbg!();

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
        // println!("{}", line_str);

        // Parse segment line format: P<tab>p_name<tab>steps
        let mut fields = line_str.split('\t');
        
        let Some((name, steps)) = fields.next().and_then(|_type| {
            let name = fields.next()?;
            let steps = fields.next()?;
            Some((name, steps))
        }) else {
            println!("Unable to parse GFA path: {}\n", line_str);
            continue;
        };

        // check if path name in map
        if !path_to_seq_id_map.contains_key(name) {
            println!("Skipping path: {} not found\n", name);
            continue;
        }

        let seq_id = path_to_seq_id_map[name] * 2;  // Paths are represented by positive strands
        let st = path_starts[seq_id];
        let end = path_starts[seq_id + 1];
        let num_matches = end - st;

        if num_matches <= 0 {
            println!("No entry for path: {}, seq_id: {}, st:{}\n", name, seq_id, st);
            continue;
        }

        let steps: Vec<&str> = steps.split(',').collect();

        println!("Processing path: {}   seq_id: {}   st: {}   length: {}", name, seq_id, st, num_matches);


        let mut processed_matches = 0;

        let mut path_pos_reader = create_reader(path_pos_path)?;
        let mut pos_line = String::new();
        let mut buffer = Vec::new();

        // Skip lines without string allocation
        for _ in 0..st {
            buffer.clear();
            path_pos_reader.read_until(b'\n', &mut buffer)?;
        }

        path_pos_reader.read_line(& mut pos_line)?;     // todo: handle empty line?
        let pos_line_str = pos_line.trim();
        let mut fields :Vec<&str> = pos_line_str.split('\t').collect();

        if fields.len() < 6 {
            eprintln!("Warning: Insufficient fields in path_pos line\n");
            continue;
        }
        let (mut curr_node_id, mut curr_offset, mut alignment_len, mut read_start, mut read_id) = (
            fields[0].parse::<usize>().unwrap(),
            fields[1].parse::<usize>().unwrap(),
            fields[3].parse::<usize>().unwrap(),
            fields[4].parse::<usize>().unwrap(),
            fields[5].parse::<usize>().unwrap(),
        );

        println!("Processing node: {} offset: {}, al_len: {}, read_id: {}, read_st: {}", curr_node_id, curr_offset, alignment_len, read_id, read_start);

        for (i, step) in steps.iter().enumerate() {
            let (seg, orient) = step.split_at(step.len() - 1);
            let seg_id = seg.parse::<usize>().unwrap();
            
            // continue walking down the path
            if seg_id != curr_node_id {
                continue;
            }

            processed_matches += 1;
            let remaining_steps = &steps[i..];

            // todo: handle reverse strand traversal
            // Start traversing starting from the current node
            let (path_str, path_len) = traverse_nodes(&remaining_steps, &mut callback, &mut get_node_len, curr_offset, alignment_len);
            let path_start = curr_offset;
            let path_end = curr_offset + alignment_len;
            if path_len < alignment_len {
                eprintln!("ERROR Path len {} << alignment_len {}. Path name: {}, Path str: {}", path_len, alignment_len, name, path_str);
                break;
            }
            if path_len < path_end {
                eprintln!("ERROR: Path len {} <= path_end {}. Path name: {}, Path str: {}", path_len, path_end, name, path_str);
                break;
            }

            println!("\nAlignment: read_id, read_st, path_str, path_len, path_st, path_end, alignment_len, path_name");
            println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n", read_id, read_start, path_str, path_len, path_start, path_end, alignment_len, name);
            if processed_matches == num_matches {
                break;
            }

            // load a new match
            pos_line.clear();
            path_pos_reader.read_line(& mut pos_line)?;     // todo: handle empty line?
            let pos_line_str = pos_line.trim();
            let mut fields :Vec<&str> = pos_line_str.split('\t').collect();
            (curr_node_id, curr_offset, alignment_len, read_start, read_id) = (
                fields[0].parse::<usize>().unwrap(),
                fields[1].parse::<usize>().unwrap(),
                fields[3].parse::<usize>().unwrap(),
                fields[4].parse::<usize>().unwrap(),
                fields[5].parse::<usize>().unwrap(),
            );
            println!("Processing node: {} offset: {}, al_len: {}, read_id: {}, read_st: {}", curr_node_id, curr_offset, alignment_len, read_id, read_start);
        }
        println!("---------------------");
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
fn read_path_name_indices(filename: &str) -> std::io::Result<HashMap<String, usize>> {
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

    // Print the contents of the path_indices HashMap
    eprintln!("Path indices map:");
    for (path_name, index) in &path_indices {
        eprintln!("  {}: {}", path_name, index);
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
    /// Path positions TSV file
    #[arg(long)]
    path_pos: Option<String>,
    /// Sequence ID starts file
    #[arg(long)]
    seq_id_starts: Option<String>,
    /// Path names
    #[arg(long)]
    path_names: Option<String>,
}

fn main() {
    let args = Args::parse();
    let gfa_file = &args.gfa;

    // Check if we're in CSV processing mode
    if let (Some(path_pos_file), Some(seq_id_starts_file)) = (&args.path_pos, &args.seq_id_starts) {
        // let path_pos_data = read_path_pos_csv(path_pos_file).unwrap();
        println!("Processing gfa file: {}", gfa_file);
        println!("Args: {}, {}", path_pos_file, seq_id_starts_file);
        let path_names_file_name = &args.path_names.unwrap();
        println!("Path names file: {}", path_names_file_name);
        // let seq_id_starts_map = read_seq_id_starts_map(seq_id_starts_file).unwrap();

        let seq_starts = read_cumulative_starts(seq_id_starts_file).unwrap();
        let path_to_seq_id_map = read_path_name_indices(path_names_file_name).unwrap();

        // Parse GFA file
        let (segment_lengths, min_id) = parse_gfa(&gfa_file).unwrap();
        let num_segments = segment_lengths.len();       // segment -> node

        let mut coverage: Vec<f64> = vec![0.0; num_segments];

        if let Err(e) = walk_gfa(&gfa_file, path_pos_file, path_to_seq_id_map, seq_starts, |node_id, len| {
            coverage[node_id - min_id] += len as f64;
        }, |node_id| segment_lengths[node_id - min_id]) {
            eprintln!("Error processing GFA file with path positions: {}", e);
            std::process::exit(1);
        }
        print!("#sample");
        for n in min_id..min_id + num_segments {
            print!("\tnode.{}", n);
        }
        println!();
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
        println!();
        return;
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
            println!(
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
        println!();
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
        println!();
    }
}
