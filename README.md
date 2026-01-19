# gafpack

A tool for working with pangenome variation graphs that:
1. Calculates node coverage from GAF alignments to GFA variation graphs
2. Generates GAF alignments by walking GFA paths and matching them with pre-computed positions

This is useful for:
- Analyzing read coverage across pangenome variation graphs
- Supporting haplotype-based genotyping workflows
- Quantifying alignment distribution across graph nodes
- Generating GAF alignments from graph path traversals

## Install

Install using:

```bash
cargo install --git https://github.com/pangenome/gafpack
```

Or build from source:

```bash
git clone https://github.com/pangenome/gafpack
cd gafpack
cargo build --release
```

## Usage

### Mode 1: Calculate Coverage from GAF Alignments

Basic usage:

```bash
gafpack --gfa graph.gfa --gaf alignments.gaf > coverage.tsv
```

The GFA file can be gzip/bgzip compressed (`.gz` or `.bgz`).

### Mode 2: Generate GAF from Path Positions (New)

Generate GAF alignments by walking GFA paths and matching them with pre-computed path positions:

```bash
gafpack --gfa graph.gfa \
        --path-pos matches.tsv \
        --seq-id-starts seq_starts.txt \
        --path-names paths.txt \
        --gaf-file-prefix output
```

This mode:
- Walks through all paths in the GFA file (both forward and reverse)
- Matches reads to path positions using pre-computed match data
- Generates GAF alignment entries for each match
- Outputs both GAF alignments and node coverage statistics

**Output files:**
- `{prefix}.gaf` - GAF format alignments with columns: read_id, read_st, path_str, path_len, path_st, path_end, match_len, path_name
- `{prefix}_coverage.csv` - Node coverage statistics (node_id, node_coverage)

## Options

### Required Arguments

- `--gfa`: Input GFA graph file (required, supports .gz/.bgz compression)

### Mode 1: GAF Coverage Analysis

- `-g, --gaf`: Input GAF alignment file (required for Mode 1)
- `-l, --len-scale`: Scale coverage by node length
- `-c, --coverage-column`: Output coverage vector as single column
- `-w, --weight-queries`: Weight coverage by query occurrences

### Mode 2: Path Position Processing

- `--path-pos`: Path positions TSV file containing match locations (tab-separated: node_id, offset, strand, mem_len, read_st, read_id)
- `--seq-id-starts`: Sequence ID starts file (cumulative match positions per path/strand)
- `--path-names`: Path names file (one path name per line, maps to sequence IDs)
- `--gaf-file-prefix`: Output file prefix for GAF and coverage files
- `-l, --len-scale`: (Optional) Scale coverage values by node length in output

## Output Formats

### Mode 1: Coverage from GAF Alignments

#### Default (tabular):

```
#sample        node.1  node.2  node.3  ...
alignments.gaf 1.5     2.0     0.5     ...
```

#### Column format (with `-c, --coverage-column`):

```
##sample: alignments.gaf
#coverage
1.5
2.0
0.5
...
```

### Mode 2: GAF Generation from Path Positions

#### GAF Output (`{prefix}.gaf`):

Tab-separated format with header:
```
read_id  read_st  path_str         path_len  path_st  path_end  match_len  path_name
1        0        >123>456>789     5000      100      300       200        chr1
2        50       <789<456<123     5000      200      400       200        chr1_reverse
```

**Columns:**
- `read_id`: Read/query identifier
- `read_st`: Start position in the read
- `path_str`: Path through the graph (e.g., `>123>456` for forward, `<789<456` for reverse)
- `path_len`: Total length of the path
- `path_st`: Start position on the path
- `path_end`: End position on the path
- `match_len`: Length of the alignment/match
- `path_name`: Name of the reference path

#### Coverage Output (`{prefix}_coverage.csv`):

```
node_id,node_coverage
1,45.50
2,23.00
3,67.25
```

## Input File Formats

### Path Positions File (`--path-pos`)

Tab-separated file with match positions:
```
node_id  offset  strand  mem_len  read_st  read_id
123      100     1       200      0        1
456      50      0       150      200      1
```

**Columns:**
- `node_id`: Graph node identifier
- `offset`: Offset within the node
- `strand`: Strand (0 or 1)
- `mem_len`: Match length
- `read_st`: Start position in the read
- `read_id`: Read identifier

### Sequence ID Starts File (`--seq-id-starts`)

File containing cumulative start positions (one per line):
```
0
5
12
```

Each line represents the cumulative number of matches up to that sequence ID (path/strand combination).

### Path Names File (`--path-names`)

One path name per line:
```
chr1
chr2
chrX
```

Line index maps to sequence ID (0-based). Each path generates two sequence IDs (forward and reverse strand).

## Examples

### Example 1: Basic Coverage Calculation

```bash
# Calculate coverage from GAF alignments
gafpack --gfa pangenome.gfa --gaf reads.gaf > coverage.tsv

# With length scaling
gafpack --gfa pangenome.gfa --gaf reads.gaf --len-scale > normalized_coverage.tsv

# With query weighting (useful for multi-mapped reads)
gafpack --gfa pangenome.gfa --gaf reads.gaf --weight-queries > weighted_coverage.tsv
```

### Example 2: Generate GAF from Path Matches

```bash
# Process pre-computed path matches and generate GAF alignments
gafpack --gfa pangenome.gfa \
        --path-pos read_matches.tsv \
        --seq-id-starts cumulative_starts.txt \
        --path-names reference_paths.txt \
        --gaf-file-prefix my_output

# This creates:
#   - my_output.gaf (alignments)
#   - my_output_coverage.csv (node coverage)
```

### Example 3: Path Processing with Length Scaling

```bash
# Generate GAF and compute length-normalized coverage
gafpack --gfa pangenome.gfa \
        --path-pos matches.tsv \
        --seq-id-starts starts.txt \
        --path-names paths.txt \
        --gaf-file-prefix results \
        --len-scale
```

### Example 4: Real-World Usage

```bash
# Process a specific genomic region with path positions
gafpack --gfa s28cc_gla_chr2_274000_284000.gfa \
        --path-pos output_path_pos.tsv \
        --seq-id-starts output_seq_id_starts.out \
        --path-names s28cc_gla_chr2_274000_284000.paths \
        --gaf-file-prefix output \
        --len-scale

# Output files created:
#   - output.gaf (alignments)
#   - output_coverage.csv (length-normalized node coverage)
```
