#if sequence is over 200kbp
# we will need to split it and run it in separate threads
#this will do it hopefully

import os
import tempfile
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from blat_main import *

# === CONFIGURATION ===
input_fasta = "/usr1/gouallin/blat/franken_plasmid_tests/franken_arg/franken_arg.fasta"                         # or full path, e.g., "/usr1/shared/pblat/blat_main"
chunk_size = 100000                          # 100kb
max_threads = 46                   # Use half of all available threads

output_dir = "franken_arg_blat_outputs_100k"
os.makedirs(output_dir, exist_ok=True)

# === CHUNKING + BLAT FUNCTION ===
def run_blat_on_chunk(args):
    species = "Enterococcus_phage"
    chunk_record, idx, original_id = args


    #try and run 114 and onwards again
    if idx <= 114:
        return None

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".fa") as tmp_fasta:
        SeqIO.write(chunk_record, tmp_fasta, "fasta")
        tmp_fasta_path = tmp_fasta.name

    output_file = os.path.join(output_dir, f"{original_id}_chunk_{idx}")
    print("output_file:", output_file)
    blat_cmd = ["python3",
                "/usr1/gouallin/blat/blat_pipeline/blat_main.py",
                tmp_fasta_path, 
                output_file, 
                "/usr1/shared/all_gtdb_seq_species_location_index.pkl",
                "/usr1/shared/TimeTree_v5_Final.nwk",
                "/usr1/shared/gtdb_split_2bit_1k/",
                species]


    print(f"Running: {' '.join(blat_cmd)}")
    try:
        subprocess.run(blat_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAT on chunk {idx}: {e}")
    finally:
        os.remove(tmp_fasta_path)
    assert(False)

# === MAIN ===
def main():
    chunk_jobs = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)
        if seq_len > chunk_size:
            for i in range(0, seq_len, chunk_size):
                chunk_seq = record.seq[i:i+chunk_size]
                chunk_record = record[:0]  # copy header
                chunk_record.seq = chunk_seq
                chunk_record.id = f"{record.id}_chunk_{i // chunk_size}"
                chunk_record.description = f"{record.description} chunk {i // chunk_size}"
                chunk_jobs.append((chunk_record, i // chunk_size, record.id))
        else:
            chunk_jobs.append((record, 0, record.id))

    print(f"Prepared {len(chunk_jobs)} chunk(s) to BLAT.")

    #this will parallelize it
    with Pool(processes=max_threads) as pool:
        pool.map(run_blat_on_chunk, chunk_jobs)

    print("All BLAT jobs finished.")

if __name__ == "__main__":
    main()