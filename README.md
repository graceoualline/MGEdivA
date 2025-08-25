# Blatdiver

A parallelized sequence alignment tool that uses BLAT and divergence to find Mobile Genetic Elements (MGEs) within sequences.

## Overview

Blatdiver is designed for optimized threading of BLAT (BLAST-Like Alignment Tool) operations. It splits large sequences into manageable chunks, runs BLAT searches in parallel across multiple database files, and applies divergence filtering and sequence analysis to identify potential MGEs.

## Installation

### Prerequisites
Please ensure you have the following tools downloaded:
- Python 3.7+
- BLAT: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/
- Kraken2: https://github.com/DerrickWood/kraken2/wiki/Manual

## Usage

### Quick Start

```bash
python3 blatdiver.py -q input.fasta -o output_directory -d /path/to/blat_database -tr TimeTree_v5_Final.nwk -i database_index.txt -k /path/to/kraken_db
```

### Parameters

#### Required Arguments
- `-q, --query`: Path to the query FASTA file
- `-o, --output`: Name of your output directory
- `-d, --database`: Path to the BLAT database directory
- `-tr, --tree`: Path to the TimeTree of Life tree file (TimeTree_v5_Final.nwk)
- `-i, --index`: Index file mapping sequences in BLAT database to species and locations
- `-k, --kraken`: Path to Kraken2 database

#### Optional Arguments
- `-t, --threads`: Number of threads to use (default: 1) # Highly recommend using more threads to speed up the program.
- `-c, --chunk`: Chunk size for sequence splitting (default: 100,000 bp)
- `-r, --remove`: Clean up intermediate files (1=True, 0=False, default: 1)
- `-s, --species`: Predefined species name (replace spaces with '_')
- `-minScore`: Minimum BLAT alignment score (default: 30)
- `-minIdentity`: Minimum percent identity threshold (default: 0)
- `-overlap_filter`: Enable overlap filtering (1=True, 0=False, default: 0)
- `-overlap_div_filter`: Enable overlap and divergence filtering (1=True, 0=False, default: 0)

### Example Commands

#### Basic run with 8 threads:
```bash
python blatdiver.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d /data/blat_db \
  -tr TimeTree_v5_Final.nwk \
  -i database_index.txt \
  -k /data/kraken2_db \
  -t 8
```

#### Run with additional filtering:
```bash
python blatdiver.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d /data/blat_db \
  -tr TimeTree_v5_Final.nwk \
  -i database_index.txt \
  -k /data/kraken2_db \
  -t 8 \
  -overlap_filter 1 \
  -overlap_div_filter 1 \
  -minIdentity 95
```

## Output Files

BlatDiver generates several output files in the specified output directory:

- `{output_name}_blat_results.tsv`: Raw BLAT alignment results
- `{output_name}_blatdiver_output.tsv`: Filtered results with divergence information
- `{output_name}_overlap.tsv`: Overlap-filtered results (if enabled)
- `{output_name}_overlap_div.tsv`: Overlap and divergence filtered results (if enabled)

## Database Requirements

### BLAT Database
- Must contain `.2bit` files and corresponding `.ooc` files
- Each `.2bit` file must have a matching `.ooc` file
- Files should be properly indexed

### Index File
- Tab-separated file mapping database sequences to species and genomic locations
- Format: `sequence_id\tspecies\tlocation_info`

### TimeTree File
- Newick format phylogenetic tree
- Should contain species present in your analysis

## Performance Tips

1. **Thread Usage**: Use multiple threads (`-t`) for faster processing, but don't exceed your CPU core count
2. **Chunk Size**: Adjust chunk size (`-c`) based on your sequence lengths and memory constraints
3. **Filtering**: Enable additional filters only if needed, as they increase processing time
4. **Memory**: Large databases may require substantial RAM; monitor system resources

## Workflow

1. **Input Processing**: Reads FASTA sequences and determines species classification
2. **Chunking**: Splits sequences larger than chunk size into manageable pieces
3. **BLAT Search**: Runs parallel BLAT searches against database files
4. **Result Combination**: Merges BLAT outputs from all database chunks
5. **Filtering**: Applies divergence, overlap, and other specified filters
6. **Output Generation**: Combines all results into final output files
7. **Cleanup**: Removes intermediate files (if enabled)

## Troubleshooting

### Common Issues

- **File Permission Errors**: Ensure read/write permissions for input/output directories
- **Memory Issues**: Reduce chunk size or number of threads
- **Missing Dependencies**: Verify all required modules are installed and accessible
- **Database Mismatch**: Ensure equal numbers of `.2bit` and `.ooc` files

### Error Messages

- `Mismatch in file counts`: Check that every BLAT file has a matching OOC file
- `Species extraction failed`: Verify Kraken2 database path and format
- `BLAT command failed`: Check BLAT installation and database integrity

## Citation

If you use BlatDiver in your research, please cite:
[Add appropriate citation information]

## License

[Add license information]

## Contributing

[Add contribution guidelines]

## Support

For questions and support, please [add contact information or issue tracker link].
