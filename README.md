# MGEdivA

A parallelized sequence alignment tool that uses BLAT and divergence to find Mobile Genetic Elements (MGEs) within sequences.

## Overview

MGEdivA is designed to find novel MGEs through optimized threading of BLAT (BLAST-Like Alignment Tool) operations, and advanced filtering. It splits large sequences into manageable chunks, runs BLAT searches in parallel across multiple database files, and applies divergence filtering and sequence analysis to identify potential MGEs. The tool is particularly effective for horizontal gene transfer detection.

## Installation
Install MGEdivA:
```bash
git clone https://github.com/graceoualline/MGEdivA.git
cd mgediva
```
Install supporting files and databases here: (PUT LINK TO LARGER FILES)
You should have the following:
- gtdb_2bil_split_2bit
- TimeTree_v5_Final.nwk
- gtdb_2bil_seq_id_species_loc_index.pkl
- kraken2_custom_db
### Required Python Packages
```bash
pip install biopython tqdm
```
### Prerequisites
Please ensure you have the following tools downloaded:
- Python 3.7+
- BLAT: https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/
  ```bash
  conda install -c bioconda blat
  ```
- Kraken2: https://github.com/DerrickWood/kraken2/wiki/Manual
  ```bash
  conda install -c bioconda kraken2
  ```
#### Kraken and Kraken Database
We provide the kraken database ```kraken2_custom_db```. If you want to recreate it locally, do the following:
```bash
mkdir -p kraken2_custom_db

kraken2-build --download-taxonomy --db kraken2_custom_db

kraken2-build --download-library archaea --db kraken2_custom_db
kraken2-build --download-library bacteria --db kraken2_custom_db
kraken2-build --download-library plasmid --db kraken2_custom_db
kraken2-build --download-library viral --db kraken2_custom_db
kraken2-build --download-library fungi --db kraken2_custom_db
kraken2-build --download-library protozoa --db kraken2_custom_db
kraken2-build --download-library nt --db kraken2_custom_db

kraken2-build --build --db kraken2_custom_db
kraken2-build --clean --db kraken2_custom_db
```
## Usage

### Quick Start

```bash
# To see all input parameters
python3 mgediva.py -h

# command line with only required arguments
python3 mgediva.py -q input.fasta -o output_directory -d /path/to/blat_database -tr TimeTree_v5_Final.nwk -i database_index.txt -k /path/to/kraken_db --threads 20

# with config
python mgediva.py --config config_example.yaml

# with config and command line. Command-line arguments take priority over config file values.
python mgediva.py --config config_example.yaml -q input.fasta -o output_directory --threads 20
```
#### Ready-to-Run Example (only on FAUST) 
```bash
# with config
python3  /usr1/gouallin/blat/blat_pipeline/mgediva.py --config /usr1/gouallin/blat/blat_pipeline/config_example.yaml

# command line only
python3 /usr1/gouallin/blat/blat_pipeline/mgediva.py \
  -q /usr1/gouallin/blat/blat_pipeline/test/acrB.fasta \
  -o mgediva_test_results_acrB \
  -d /usr1/shared/gtdb_2bil_split_2bit/ \
  -tr /usr1/shared/TimeTree_v5_Final.nwk \
  -i /usr1/shared/gtdb_2bil_seq_id_species_loc_index.pkl \
  -k /usr1/shared/kraken2_custom_db/ \
  -t 20 \
  --no-remove \
  -minIdentity 90 \
  --overlap_filter \
  --overlap_div_filter
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
- `-t, --threads`: Number of threads to use (default: 1) - **Highly recommend using more threads to speed up the program.**
- `-c, --chunk`: Chunk size for sequence splitting (default: 100,000 bp)
- `-s, --species`: Predefined species name (replace spaces with '_'). For multifasta files, applies to all sequences.
- `-minScore`: Minimum BLAT alignment score (default: 30)
- `-minIdentity`: Minimum percent identity threshold (default: 0)
- `--remove, --no-remove`: Clean up temporary files (default: True). Use --remove to enable, --no-remove to disable.
- `--overlap-filter, --no-overlap-filter`: Enable overlap filter (default: False). Use --overlap-filter to enable, --no-overlap-filter to disable.
- `--overlap-div-filter, --no-overlap-div-filter`: Enable overlap and divergence filter (default: False). Use --overlap-div-filter to enable, --no-overlap-div-filter to disable.

#### Config File
There is an option for a config file, as the majoriy of parameters are likely reused in each run. An example config file is available in `config_example.yaml`. Here is an example config file:
```
# Example config, for MGEdivA.
# Required arguments:
query: test/acrB.fasta
output: mgediva_test_results_acrB
database: gtdb_2bil_split_2bit/
index: gtdb_2bil_seq_id_species_loc_index.pkl
kraken: kraken2_custom_db/ 
tree: TimeTree_v5_Final.nwk 

# Optional Arguments
chunk: 100000
threads: 46
remove: true
species: null
minScore: 30
minIdentity: 95 
overlap_filter: false
overlap_div_filter: false
```

### Example Commands

#### Basic Command with GTDB-Blat Database and 20 threads
```bash
python mgediva.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d gtdb_2bil_split_2bit \
  -tr TimeTree_v5_Final.nwk \
  -i gtdb_2bil_seq_id_species_loc_index.pkl \
  -k kraken2_custom_db \
  -t 20 \
```

#### Run with additional filtering:
```bash
python mgediva.py \
  -q my_sequences.fasta \
  -o results_dir \
  -d gtdb_2bil_split_2bit \
  -tr TimeTree_v5_Final.nwk \
  -i  gtdb_2bil_seq_id_species_loc_index.pkl \
  -k kraken2_custom_db \
  -t 20 \
  --overlap_filter \
  --overlap_div_filter \
  -minIdentity 95
```
## Filters
Filters are described in further detail, and their processes are illustrated in our paper (add cite).
### Divergence Filtering
- Divergence filtering is the main method for detecting MGEs and produces the file ```{output_name}_mgediva_output.tsv```.
- By examining the species of the query genome and the genome it aligned to, we use the TimeTree of Life to calculate the divergence time between the two species. If the species diverged over 1 million years ago (```divergence >= 1 MYA```), the alignment is retained.
- This filter is effective at identifying horizontal gene transfer events because MGEs transferred between distantly related species will show high sequence similarity despite ancient species divergence.
- For detailed information on how this filter detects MGEs, please refer to our paper: (citation tba).
### Additional filtering
These filters further refine the alignments that were identified as divergently distant in ```{output_name}_mgediva_output.tsv```:
#### Overlap Filtering
- Enabling ```--overlap_filter``` produces the file ```{output_name}_overlap.tsv```
- This filter processes the results from ```{output_name}_mgediva_output.tsv```
- It identifies two alignments that overlap spatially on the query sequence but originate from different species
- While effective at reducing false positives, it also reduces mgediva's sensitivity for detecting true MGEs (see paper for details)
#### Overlap Divergence Filtering
- Enabling ```--overlap_div_filter``` produces the file ```{output_name}_overlap_div.tsv```
- This filter processes the results from {output_name}_mgediva_output.tsv
-This filter identifies overlapping alignments where the source species are divergently distant (>= 1 MYA) 
- It is the most stringent filter that effectively reduces false positives but may decrease sensitivity for detecting true MGEs and requires longer processing time.
- Recommended for high-confidence MGE detection when processing time is not a constraint

## Output Files
mgediva generates several output files in the specified output directory:

- `{output_name}_blat_results.tsv`: Raw BLAT alignment results
- `{output_name}_mgediva_output.tsv`: Filtered results with divergence information
- `{output_name}_overlap.tsv`: Overlap-filtered results (if enabled)
- `{output_name}_overlap_div.tsv`: Overlap and divergence filtered results (if enabled)

### Resume Functionality
Important: The program is designed to resume from interruptions by checking for existing files. If a run is stopped prematurely, it will restart from where it left off. Avoid creating files with names that could overlap with mgediva's output to prevent conflicts.

## Database Setup

### Pre-build GTDB-Blat Database
We provide a ready-to-use BLAT-compatible database created from the Genome Taxonomy Database (GTDB):
- **Database:** ```gtdb_2bil_split_2bit/```
- **Index:** ```gtdb_2bil_seq_id_species_loc_index.pkl```

### Creating a Custom Blat database
-  If you wish to create your own database, follow these steps:
#### Step 1: Build the Blat database
The database must be in 2bit format with a maximum of ~2 billion base pairs per file. We recommend a size of 2 billion base pairs.
```
python3 make_blat_db.py [fasta_file] [output_name] [size_in_bil_bp]
```
For example:
```
python3 make_blat_db.py gtdb.fa blat_gtdb_db 2
```
This script will:
- Split the multifasta file into chunks of specified size
- Convert split files to 2bit format
- Generate .ooc files for all 2bit files
- Ensure each .2bit file has a matching .ooc file
#### Step 2: Extract Species Information
Classify sequences using Kraken2:
```
python3 extract_species_from_kraken.py [fasta file] [kraken_db]
```
##### Outputs:
- ```{fasta_file}_kraken_output.txt```: Kraken's detailed classification output
- ```{fasta_file}_seq_species.txt```: Tab-separated file iwth sequence ID and species:
```
seq_id1  species1
seq_id2  species2
...  ...
```
You can manually edit species assignments in the ```{fasta_file}_seq_species.txt``` file if needed.
#### Step 3: Build the Database Index
Create a hash table index linking sequence IDs to their location and species:
```
python3 build_database_index.py [blat_db_2bit_dir] [output_index_name.pkl] [{fasta_file}_seq_species.txt]
```
This creates an index where ```index[seq_id] = (species, location)```.

## Required Files Summary
1. **BLAT Database:** Directory containing .2bit and .ooc files
2. **Index File:** .pkl file mapping sequence IDs to species and locations
3. **TimeTree File:** Newick format phylogenetic tree with relevant species
4. **Kraken Database:** For taxonomic classification

## Performance Tips
1. Use Multiple Threads: Set -t to utilize available CPU cores (recommended: 10-20 threads)
2. Optimal Chunk Size: Default 100kb works well; adjust based on sequence lengths. Increasing chunk size will significantly down the program.
3. Enable Filtering Judiciously: Enabling overlap and overlap divergence filters will decrease false positives, but will also cause mgediva ot miss more MGEs. Overlap divergence takes a long time to run.
4. Monitor Resources: Large databases require substantial RAM
5. Resume Feature: Take advantage of the resume capability for long runs

## Workflow
**Input Processing:** Reads FASTA sequences and (if not specified by the user) determines species classification via Kraken2
**Chunking:** Splits sequences larger than chunk size into manageable pieces
**Parallel BLAT Search:** Runs BLAT searches against database files using multiple threads
**Filtering Pipeline:** Applies divergence other specified filters
**Cleanup:** Removes intermediate files (if enabled)

## Citation

If you use mgediva in your research, please cite:
[Add later]

## Support

For questions and support, please submit an issue ticket on the GitHub repository.
