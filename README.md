# Blatdiver

A parallelized sequence alignment tool that uses BLAT and divergence to find Mobile Genetic Elements (MGEs) within sequences.

## Overview

Blatdiver is designed for optimized threading of BLAT (BLAST-Like Alignment Tool) operations. It splits large sequences into manageable chunks, runs BLAT searches in parallel across multiple database files, and applies divergence filtering and sequence analysis to identify potential MGEs.

## Installation
```
git clone https://github.com/graceoualline/blatdiver.git
```
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
#### Basic ready-to-run commmand for Dr.Yu 
```bash
python /usr1/gouallin/blat/blat_pipeline/blatdiver.py \
  -q /usr1/gouallin/blat/blat_pipeline/test/acrB.fasta \
  -o blatdiver_test_results_acrB \
  -d /usr1/shared/gtdb_2bil_split_2bit/ \
  -tr /usr1/shared/TimeTree_v5_Final.nwk \
  -i /usr1/shared/gtdb_2bil_seq_id_species_loc_index.pkl \
  -k /usr1/shared/kraken2_custom_db/ \
  -t 20 \
  -r 0 \
  -minIdentity 90 \
  -overlap_filter 1 \
  -overlap_div_filter 1
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
- `-s, --species`: If you want to define the species of your query, otherwise we will find it for you with Kraken. (replace spaces with '_') (For multifasta files, we will apply this species to all of the sequences in the file. We do not have the ability to manually assign different species to each individual sequence.)
- `-minScore`: Minimum BLAT alignment score (default: 30)
- `-minIdentity`: Minimum percent identity threshold (default: 0)
- `-overlap_filter`: Enable overlap filtering (1=True, 0=False, default: 0)
- `-overlap_div_filter`: Enable overlap and divergence filtering (1=True, 0=False, default: 0)

### Example Commands

#### Basic Command
```bash
python blatdiver.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d /usr1/shared/gtdb_2bil_split_2bit/ \
  -tr /usr1/shared/TimeTree_v5_Final.nwk \
  -i /usr1/shared/gtdb_2bil_seq_id_species_loc_index.pkl \
  -k /usr1/shared/kraken2_custom_db/ \
  -t 20 \
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
  -t 20 \
  -overlap_filter 1 \
  -overlap_div_filter 1 \
  -minIdentity 95
```

## Output Files

blatdiver generates several output files in the specified output directory:

- `{output_name}_blat_results.tsv`: Raw BLAT alignment results
- `{output_name}_blatdiver_output.tsv`: Filtered results with divergence information
- `{output_name}_overlap.tsv`: Overlap-filtered results (if enabled)
- `{output_name}_overlap_div.tsv`: Overlap and divergence filtered results (if enabled)

*Note: The program is designed such that if a run is stopped prematurely, it will be able to start again at the point where it left off. It does this by checking for what files already exist. Therefore, try to avoid creating any files with names that could overlap with blatdiver's output.

## Database Requirements

### BLAT Database
- We have created a blat compatible database of the GTDB, named "gtdb_2bil_split_2bit".
### Creating the Blat database
-  If you wish to create your own database, it must fulfill the following requirements:
-   Each `.2bit` file must have a matching `.ooc` file
-   Files should have a proper index made that will contain each sequence's location within the database, and its species (instructions below).

#### Building the database
- The blat database must be in 2bit form, and can have a maximum of around 2 billion base pairs in each 2bit file. Put all of the sequences that you want in the database into the same multifasta file. The program make_blat_db.py will do the following:
-   split the multifasta file such that ~2 billion (or specified amount) bp are in each split file.
-   convert split files to 2bit format
-   make ooc files for all of the 2bit files.
```
python3 make_blat_db.py [fasta file] [output name] [size of split in bil. bp (recommend 2)]
```
For example:
```
python3 make_blat_db.py gtdb.fa blat_gtdb_db 2
```

### Index File
- The index file is a hash table that contains all sequence ids in the database. It allows us to quickly reference sequence's location and species during the pipeline.
- The index file created for the blat-gtdb database is named "gtdb_2bil_seq_id_species_loc_index.pkl".

#### Creating an index file for a custom database
- First, you will need to run kraken on all of your sequences, so that we are able to classify their species.
```
python3 extract_species_from_kraken [fasta file] [kraken_db]
```
- The output {fasta_file}_kraken_output.txt, which is the output that kraken generates showing the tree splits it made to characterize your sequence
- The output {fasta_file}_seq_species.txt is a tab separated file like so:
```
seq_id1  species1
seq_id2  species2
```
- If you wish to manually define some of the sequences in your database, you can edit the species next to that sequence in {fasta_file}_seq_species.txt
- The script build_database_index.py will intake a blat database of 2bit files and {fasta_file}_seq_species.txt, and output an index linking sequence ids to their location and species. It will do the following:
-   get the locations of all of the sequences within the 2bit files
-   creates a hash table index where index[seq_id] = (species, location)
```
python3 build_database_index.py [blat_db_2bit_dir] [output_index_name.pkl] [{fasta_file}_seq_species.txt]
```

### TimeTree File
- Newick format phylogenetic tree
- Should contain species present in your analysis

### Kraken Database
- f

## Troubleshooting
- 

## Common Issues

- 

## Error Messages

- 

## Citation

If you use blatdiver in your research, please cite:
[Add later]

## Support

For questions and support, please submit an issue ticket.
