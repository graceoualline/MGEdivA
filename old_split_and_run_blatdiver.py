#if sequence is over 200kbp
# we will need to split it and run it in separate threads
#this will do it hopefully

import os
import tempfile
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from blat_main import *
from datetime import datetime
from extract_species_from_kraken import *

# === CONFIGURATION ===                  
#max_threads = 46                   # Use half of all available threads
'''
# full franken_arg
input_fasta = "/usr1/gouallin/blat/franken_plasmid_tests/franken_arg/franken_arg.fasta"
chunk_size = 100000 
output_dir = "franken_arg_blat_outputs_100k"
species = "Enterococcus_phage"
database = "/usr1/shared/gtdb_split_2bit_1k/"
index = "/usr1/shared/all_gtdb_seq_species_location_index.pkl"

#normal database split 1000
input_fasta = "/usr1/gouallin/blat/franken_plasmid_tests/test_arg_results/one_lucky_arg.fasta"
chunk_size = 1000
output_dir = "one_lucky_arg_blat_outputs_1000"
species = "Salmonella_enterica"
database = "/usr1/shared/gtdb_split_2bit_1k/"
index = "/usr1/shared/all_gtdb_seq_species_location_index.pkl"

#smart_split database
input_fasta = "/usr1/gouallin/blat/franken_plasmid_tests/test_arg_results/one_lucky_arg.fasta"
chunk_size = 1000
output_dir = "one_lucky_arg_blat_outputs_1000_smart_split"
species = "Salmonella_enterica"
database = "/usr1/gouallin/blat/smart_split_db/gtdb_smart_split_2bit"
index = "/usr1/gouallin/blat/smart_split_db/all_gtdb_seq_species_location_index_smart_split.pkl"

#franken_100_plasmid
input_fasta = "/usr1/gouallin/blat/franken_plasmid_tests/franken_100_plasmid/franken_100_seq/franken_100_plasmid/franken_100_plasmid.fasta"
chunk_size = 100000
output_dir = "franken_100_plasmid_blatdiver"
species = "Enterococcus_phage"
database = "/usr1/shared/gtdb_split_2bit_1k/"
index = "/usr1/shared/all_gtdb_seq_species_location_index.pkl"
'''


def run_all_filters(args):
    chunk_record, idx, original_id, species = args

    #get the name and path for the output file
    output_file = os.path.join(output_dir, f"chunk_{idx}_{original_id}_blatdiver_output.tsv")
    #run the output on all of the filters

    #no gaps
    #python3 remove_large_gaps.py <blatdiver_output.tsv> <output file name>
    no_gaps = os.path.join(output_dir, f"chunk_{idx}_{original_id}_no_gaps.tsv")
    
    if os.path.exists(no_gaps):
        print(f"Skipping {no_gaps} , already exists.")
    else:
        no_gaps_cmd = ["python3",
                    "/usr1/gouallin/blat/blat_pipeline/remove_large_gaps.py", 
                    output_file,
                    no_gaps, 
                    ]
        try:
            #run the command
            print("running", no_gaps_cmd)
            subprocess.run(no_gaps_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running no gaps on chunk {idx}: {e}")
    
    #overlap
    #Usage: python3 find_overlap.py <file of blatdiver output> <output file name>
    overlap = os.path.join(output_dir, f"chunk_{idx}_{original_id}_overlap.tsv")
    if os.path.exists(overlap):
        print(f"Skipping {overlap} , already exists.")
    else:
        overlap_cmd = ["python3",
                    "/usr1/gouallin/blat/blat_pipeline/find_overlap.py",
                    output_file, 
                    overlap, 
                    database,
                    index
                    ]
        try:
            #run the command
            print("running", overlap_cmd)
            subprocess.run(overlap_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running overlap on chunk {idx}: {e}")

    #overlap and gap
    overlap_gap = os.path.join(output_dir, f"chunk_{idx}_{original_id}_overlap_and_gap.tsv")
    if os.path.exists(overlap_gap):
        print(f"Skipping {overlap_gap} , already exists.")
    else:
        overlap_gap_cmd = ["python3",
                    "/usr1/gouallin/blat/blat_pipeline/find_overlap.py",
                    no_gaps, 
                    overlap_gap,
                    "/usr1/shared/gtdb_split_2bit_1k/",
                    "/usr1/shared/all_gtdb_seq_species_location_index.pkl"
                    ]
        try:
            #run the command
            print("running", overlap_cmd)
            subprocess.run(overlap_gap_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running overlap and gap on chunk {idx}: {e}")

    #overlap and gap and divergence
    #python3 find_overlap_and_div.py <file of blatdiver output> <output file name> <tree> <blat_db> <kraken>
    
    overlap_gap_div = os.path.join(output_dir, f"chunk_{idx}_{original_id}_no_gap_overlap_div.tsv")
    if os.path.exists(overlap_gap_div):
        print(f"Skipping {overlap_gap_div} , already exists.")
    else:
        overlap_gap_div_cmd = ["python3",
                    "/usr1/gouallin/blat/blat_pipeline/find_overlap_and_div.py",
                    no_gaps, 
                    overlap_gap_div, 
                    "/usr1/shared/TimeTree_v5_Final.nwk",
                    "/usr1/shared/gtdb_split_2bit_1k/",
                    "/usr1/shared/all_gtdb_seq_species_location_index.pkl"
                    ]
        try:
            #run the command
            print("running", overlap_gap_div_cmd)
            subprocess.run(overlap_gap_div_cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error running overlap, gap, div on chunk {idx}: {e}")
    '''
    #overlap and divergence
    #this was already done for us in filter_blat_main
    overlap_div = os.path.join(output_dir, f"chunk_{idx}_{original_id}_overlap_div.tsv")
    overlap_div_cmd = ["python3",
                "find_overlap_and_div.py",
                output_file, 
                overlap_div, 
                "/usr1/shared/TimeTree_v5_Final.nwk",
                "/usr1/shared/gtdb_split_2bit_1k/",
                "/usr1/shared/all_gtdb_seq_species_location_index.pkl"
                ]
    try:
        #run the command
        subprocess.run(overlap_div_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running overlap and div on chunk {idx}: {e}")
    '''

# === CHUNKING + BLAT FUNCTION ===
def run_blat_on_chunk(args):
    chunk_record, idx, original_id, species = args
    print(f"BLATTING CHUNK {idx}")
    
    # make a temporary file that isolates the fasta on its own
    with tempfile.NamedTemporaryFile(prefix=f"chunk_{idx}_{original_id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
        SeqIO.write(chunk_record, tmp_fasta, "fasta")
        tmp_fasta_path = tmp_fasta.name

    #get the name and path for the output file
    output_file = os.path.join(output_dir, f"chunk_{idx}_{original_id}")
    print("output_file:", output_file)
    #then set the command on that fasta file
    blat_cmd = ["python3",
                "/usr1/gouallin/blat/blat_pipeline/blat_main.py",
                tmp_fasta_path, 
                output_file, 
                index,
                "/usr1/shared/TimeTree_v5_Final.nwk",
                database,
                species]
    #print whatever command you are running
    print(f"Running: {' '.join(blat_cmd)}")
    try:
        #run the command
        subprocess.run(blat_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAT on chunk {idx}: {e}")
    finally:
        #remove the temp fasta file
        os.remove(tmp_fasta_path)
    
    run_all_filters(args)

def combine_all_results(output_dir):
    #once all jobs are finished, combine all of the filters into one place:
    #python3 franken_plasmid_code/combine_blatdiver_output.py franken_arg_blat_outputs_100k 100000 _no_gaps_and_overlapping.tsv all_franken_arg_blatdiver_no_gaps_overlap.tsv
    parent_dir = os.path.dirname(output_dir)
    output_name = os.path.basename(output_dir)
    #python3 combine_blatdiver_output.py <blatdiver_output_directory> <chunk_size> <which files to combine ex. *blatdiver_output.tsv> <output_file>
    input_output = [("output_overlap_div_filtered.tsv", f"all_{output_name}_overlap_div.tsv"), 
                    ("blatdiver_output.tsv", f"all_{output_name}_blatdiver_output.tsv"),
                    ("no_gap_overlap_div.tsv", f"all_{output_name}_no_gaps_overlap_div.tsv"),
                    ("no_gaps.tsv", f"all_{output_name}_no_gaps.tsv"), 
                    ("overlap_and_gap.tsv", f"all_{output_name}_overlap_and_gap.tsv"),
                    ("overlap.tsv", f"all_{output_name}_overlap.tsv")]
    
    ##input_output = [ ("blatdiver_output.tsv", f"all_{output_dir}_blatdiver_output.tsv"), 
     #               ("no_gaps.tsv", f"all_{output_dir}_no_gaps.tsv")]
    for into, output in input_output:
        new_path = os.path.join(parent_dir, output)
        cmd = ["python3",
                "/usr1/gouallin/blat/franken_plasmid_tests/franken_plasmid_code/combine_blatdiver_output.py",
                output_dir,
                str(chunk_size),
                into,
                new_path]
        #print whatever command you are running
        print(f"Running: {' '.join(cmd)}")
        try:
            #run the command
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error combining results: {e}")
    
# === MAIN ===
if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("Usage: python3 split_and_run_blatdiver.py <input_fasta> <chunk_size> <output_dir_name> <species> <database> <index> <max_threads> <kraken_db>")
        sys.exit(1) 

    input_fasta = sys.argv[1]
    chunk_size = int(sys.argv[2]) 
    output_dir = sys.argv[3]
    species = sys.argv[4]
    database = sys.argv[5]
    index = sys.argv[6]
    max_threads = int(sys.argv[7])                   # Use half of all available threads
    kraken_db = sys.argv[8]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Start time:", datetime.now())
    
    chunk_jobs = []

    # will go through each sequence in the fasta
    for record in SeqIO.parse(input_fasta, "fasta"):
        #get the length of the sequence
        seq_len = len(record.seq)

        if species == "unk":
            try:
                with tempfile.NamedTemporaryFile(prefix=f"{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                    SeqIO.write(record, tmp_fasta, "fasta")
                    tmp_fasta_path = tmp_fasta.name
                species = get_q_species(tmp_fasta_path, kraken_db)
                os.remove(tmp_fasta_path)
            except:
                print("Issue extracting species for", record[:0])
                species = "unclassified"
            #if the sequence length is greater than the chunk size we want
        if seq_len > chunk_size:
            for i in range(0, seq_len, chunk_size):
                #then we chunk it by size
                chunk_seq = record.seq[i:i+chunk_size]
                chunk_record = record[:0]  # copy header
                #set the sequence to be this chunk
                chunk_record.seq = chunk_seq
                #set the id to be the chunk number
                #if it is the first chunk, i = 0, chunk = 0, if the first chunk, i = 100k so chunk = 1
                chunk_record.id = f"{record.id}_chunk_{i // chunk_size}"
                chunk_record.description = f"{record.description} chunk {i // chunk_size}"
                #record the actual sequence and its record, its chunk number, and its id
                chunk_jobs.append((chunk_record, i // chunk_size, record.id, species))
        else:
            #if it is smaller than the chunk size, just add it
            chunk_jobs.append((record, 0, record.id, species))

    print(f"Prepared {len(chunk_jobs)} chunk(s) to BLAT.")

    #this will parallelize it
    #run the function run_blat_on_chunk on all different threads
    print(f"running on {max_threads} threads")
    with Pool(processes=max_threads) as pool:
        pool.map(run_blat_on_chunk, chunk_jobs)

    combine_all_results(output_dir)
    print("All BLAT jobs finished.")
    print("Endtime time:", datetime.now())