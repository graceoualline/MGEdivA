# blatdiver for optimized threading
# will run all blat db on different threads
# then run filters on all results on different threads
# combine all results, filters and blat output into one file
# remove individual chunk files (if sequence over 100kbp)
# or if there are differerent sequences


import os
import tempfile
import subprocess
import multiprocessing as mp
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from blat_main import *
from datetime import datetime
from functools import partial
from extract_species_from_kraken import *
import argparse
from tqdm import tqdm
from dataclasses import dataclass
from pathlib import Path
from remove_large_gaps import *
from find_overlap import *
from combine_blatdiver_output import *
from find_overlap_and_div import *

@dataclass
class Config:
    input_fasta: str #query fasta
    chunk_size: int # how we want to split up the query to run in blat
    output_dir: str # output dir name
    database: Path # blat database
    index: str #index that maps sequences in the database to their species and locations
    max_threads: int # number of threads the user wants to use
    kraken_db: Path # path to the kraken database
    tree: str # the time tree of life tree
    remove: bool # is false when I am debugging and dont want files to be deleted
    species: str # user defined species of all of the sequences given in a query
    minScore: int # Is the minimum alignment score threshold, sets blat's parameter -minScore, default 30
    minIdentity: int # Is the threshold of the minimum percent identity a hit must have. Default: 95. Percent identity = ( match / Q_end - Q_start )*100
    overlap_filter: bool # Is false when you do not want to also do an overlap filter
    overlap_div_filter: bool # is false when you do not want to do an overlap and divergence filter

def combine_and_cleanup_psl_files(args):
    header_lines = 5
    output_chunk_dir, combined_path, remove = args

    if os.path.exists(combined_path):
        print(f"Skipping {combined_path}, already exists.")
    else:
        input_files = glob.glob(os.path.join(output_chunk_dir, "*.psl")) # test_blatdiver/chunk_0_emrk/*.psl

        with open(combined_path, "w") as outfile:
            outfile.write("match\tmismatch\trep. match\tN's\tQ gap count\tQ gap bases\tT gap count\tT gap bases\tstrand\tQ name\tQ size\tQ start\t Q end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\n")
            for fname in input_files:
                with open(fname) as infile:
                    lines = infile.readlines()
                    # Now find where the alignment lines begin
                    outfile.writelines(lines[header_lines:])
                if remove: os.remove(fname)
        if remove: os.rmdir(output_chunk_dir)

        #print(f"Combined {len(input_files)} files into {combined_path}.")

def run_div_filter(job):
    all_blat_raw, output_file, species, index, tree, tmp_fasta_path, database, minIdentity = job
    index = load_hash_table(index)
    print("Running", all_blat_raw)
    try:
        #run the first round of filtering on the raw blat data
        filter_blat(all_blat_raw, output_file, species, index, tree, tmp_fasta_path, database, minIdentity)
    except subprocess.CalledProcessError as e:
        print(f"Error running filter_blat on chunk {output_file}:\n{e.stderr.decode().strip()}")
    return 1

def run_specified_filter(job):
    """Dispatch the correct function based on job type."""
    job_type, args = job
    try:
        if job_type == "overlap":
            infile, outfile, database, index = args
            index = load_hash_table(index)
            print("overlap", outfile)
            rows = compress(infile)
            find_overlap(rows, outfile, database, index)
        elif job_type == "overlap_div":
            to_filter_file, overlap_div_file, tree, database, index = args
            index = load_hash_table(index)
            print("overlap_div", overlap_div_file)
            rows = compress(to_filter_file)
            print("compressed", len(rows))
            print(to_filter_file, overlap_div_file)
            find_overlap_and_div(rows, overlap_div_file, tree, database, index)
    except subprocess.CalledProcessError as e:
        print(f"[{job_type}] Error on {args[1] if len(args) > 1 else 'unknown'}:\n{e.stderr.decode().strip()}")
    return 1

def run_all_filters(chunk_jobs, c: Config):
    blat_combine_jobs = []
    jobs = []
    for args in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = args
        output_chunk_dir = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}")
        
        #make one big file that combines all of the blat outputs
        all_blat_raw = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if not(os.path.exists(all_blat_raw)): blat_combine_jobs.append((output_chunk_dir, all_blat_raw, c.remove))
        
        #then filter all blat results and add divergence or ani
        output_file = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blatdiver_output.tsv")
        jobs.append((chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c))

    output_name = os.path.basename(c.output_dir)
    blat_final_name = os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")
    overlap_div_final_name = os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")
    blatdiver_final_name = os.path.join(c.output_dir, f"{output_name}_blatdiver_output.tsv")
    overlap_final_name = os.path.join(c.output_dir, f"{output_name}_overlap.tsv")
    
    #first, combine all blat files
    if os.path.exists(blat_final_name): print(f"Skipping {blat_final_name}, already exists.")
    else:
        if len(blat_combine_jobs) > 0:
            with Pool(processes=c.max_threads) as pool:
                #Progress bar
                with tqdm(total=len(blat_combine_jobs), desc="Combining blat data") as pbar:
                    for _ in pool.imap_unordered(combine_and_cleanup_psl_files, blat_combine_jobs):
                        pbar.update()
        print("Finished combining blat data.")
    
    #run divergence filter
    div_jobs = []
    if os.path.exists(blatdiver_final_name): print(f"Skipping {blatdiver_final_name}, already exists.")
    else:
        for j in jobs:
            chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = j
            if os.path.exists(output_file): print(f"Skipping {output_file}, already exists.")
            else: div_jobs.append((all_blat_raw, output_file, species, c.index, c.tree, tmp_fasta_path, c.database, c.minIdentity))
        if len(div_jobs) > 0:
            #print whatever command you are running
            with Pool(processes=c.max_threads) as pool:
                #Progress bar
                with tqdm(total=len(div_jobs), desc="Running chunks through first divergence filter") as pbar:
                    for _ in pool.imap_unordered(run_div_filter, div_jobs):
                        pbar.update()
        print("Finished running divergence filter.")

    #run overlap, overlap_div filters
    make_overlap = c.overlap_filter
    make_overlap_div = c.overlap_div_filter
    if make_overlap and os.path.exists(overlap_final_name): 
        print(f"Skipping {overlap_final_name}, already exists.")
        make_overlap = False
    if make_overlap_div and os.path.exists(overlap_div_final_name): 
        print(f"Skipping {overlap_div_final_name}, already exists.")
        make_overlap_div = False
    filter_jobs = []
    for j in jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, blatdiver_output_file, c = j

        #overlap filter
        if make_overlap:
            overlap = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_overlap.tsv")
            if os.path.exists(overlap): print(f"Skipping {overlap}, already exists.")
            else: filter_jobs.append(("overlap", (blatdiver_output_file, overlap, c.database, c.index)))

        # overlap div filter
        if make_overlap_div:
            overlap_div = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_overlap_div.tsv")
            if os.path.exists(overlap_div): print(f"Skipping {overlap_div}, already exists.")
            else: filter_jobs.append(("overlap_div", (blatdiver_output_file, overlap_div, c.tree, c.database, c.index)))
    
    #run this first round of filtering
    if len(filter_jobs) > 0:
        with Pool(processes=c.max_threads) as pool:
            #Progress bars
            with tqdm(total=len(filter_jobs), desc="Running chunks through overlap, overlap div, filters") as pbar:
                for _ in pool.imap_unordered(run_specified_filter, filter_jobs):
                    pbar.update()
    print("All filters complete")
    
def run_blat(args):
    blat_dir, blat_file, ooc_file, query, output_path, minScore = args
    # Construct the command
    
    command = (
        f"blat {blat_dir / blat_file} {query} "
        f"-ooc={blat_dir / ooc_file} -tileSize=11 -minScore={minScore} "
        f"{output_path} -q=dna -t=dna"
    )
    try:
        #run the command
        #print(command)
        subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAT on {output_path}:\n{e.stderr.decode().strip()}")
    return 1

def run_blat_on_chunk(chunk_jobs, blat_files, ooc_files, c: Config):
    num_blat_db = len(blat_files)

    jobs = []
    output_name = os.path.basename(c.output_dir)
    final_blat_output = os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")
    # for each chunk, make all of the blat search jobs
    if os.path.exists(final_blat_output): 
        print(f"Skipping BLAT, {final_blat_output} already exists.")
        return
    for chunk in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = chunk
    
        # make a temporary file that isolates the fasta on its own
        #get the name and path for the output file
        all_blat_output = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if os.path.exists(all_blat_output): print(f"Skipping BLAT on chunk_{idx}_{original_id}, chunk_{idx}_{original_id}_blat_output.tsv already exists.")
        else:
            output_chunk = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}")
            if not os.path.exists(output_chunk): os.makedirs(output_chunk)
            
            for i in range(num_blat_db):
                output_path = f"chunk_{idx}_{original_id}_part_{i}.psl"
                output_path = os.path.join(output_chunk, output_path)
                if os.path.exists(output_path): print(f"{output_path} skipped, already exists")
                else:
                    jobs.append((c.database, blat_files[i], ooc_files[i], tmp_fasta_path, output_path, c.minScore))
    if len(jobs) > 0:
        #print whatever command you are running
        print("THREADS", c.max_threads, type(c.max_threads))
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(jobs), desc="Running Sequence(s) through BLAT") as pbar:
                for _ in pool.imap_unordered(run_blat, jobs):
                    pbar.update()

def run_combine_results(job):
    output_dir, chunk_size, to_combine, output_path, remove = job
    try: 
        to_remove = adjust_and_merge_tsvs(output_dir, chunk_size, output_path, to_combine)
        if remove: 
            for f in to_remove: os.remove(f)
    except subprocess.CalledProcessError as e:
        print(f"Error combining results {output_path}: {e}")
    print("Done")
    return 1

def combine_all_results(c: Config):
    #once all jobs are finished, combine all of the filters into one place:
    parent_dir = os.path.dirname(c.output_dir)
    output_name = os.path.basename(c.output_dir)
    #python3 combine_blatdiver_output.py <blatdiver_output_directory> <chunk_size> <which files to combine ex. *blatdiver_output.tsv> <output_file>
    input_output = [("blat_output.tsv", os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")),
                    ("blatdiver_output.tsv", os.path.join(c.output_dir, f"{output_name}_blatdiver_output.tsv"))]
    if c.overlap_div_filter:
        input_output.append(("overlap_div.tsv", os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")))
    if c.overlap_filter:
        input_output.append(("overlap.tsv", os.path.join(c.output_dir, f"{output_name}_overlap.tsv")))

    jobs = []
    for to_combine, output in input_output:
        if not os.path.exists(output): jobs.append((c.output_dir, c.chunk_size, to_combine, output, c.remove))
    if len(jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(jobs), desc="Combining all results") as pbar:
                for _ in pool.imap_unordered(run_combine_results, jobs):
                    pbar.update()

def process_single_record(record, config):
    chunk_jobs = []
    tmp_fastas = []
    
    # Get the length of the sequence
    seq_len = len(record.seq)
    
    # Get the species of the sequence
    if config.species is None:
        try:
            with tempfile.NamedTemporaryFile(prefix=f"{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                SeqIO.write(record, tmp_fasta, "fasta")
                tmp_fasta_path = tmp_fasta.name
            species = get_q_species(tmp_fasta_path, config.kraken_db)
            os.remove(tmp_fasta_path)
        except Exception as e:
            print(f"Issue extracting species for {record.id}: {e}")
            species = "unclassified"
    else:
        species = config.species
    
    print(f"Species of {record.id}: {species}")
    
    # If the sequence length is greater than the chunk size we want
    if seq_len > config.chunk_size:
        for i in range(0, seq_len, config.chunk_size):
            # Then we chunk it by size
            chunk_seq = record.seq[i:i+config.chunk_size]
            chunk_record = record[:0]  # copy header
            # Set the sequence to be this chunk
            chunk_record.seq = chunk_seq
            # Set the id to be the chunk number
            idx = i // config.chunk_size
            chunk_record.id = f"{record.id}_chunk_{idx}"
            chunk_record.description = f"{record.description} chunk {idx}"
            
            # Get a tmp file of the sequence
            with tempfile.NamedTemporaryFile(prefix=f"chunk_{idx}_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                SeqIO.write(chunk_record, tmp_fasta, "fasta")
                tmp_fasta_path = tmp_fasta.name
                tmp_fastas.append(tmp_fasta_path)
            
            # Record the actual sequence and its record, its chunk number, and its id
            chunk_jobs.append((chunk_record, idx, record.id, species, tmp_fasta_path))
    else:
        # Get a tmp file of the sequence
        with tempfile.NamedTemporaryFile(prefix=f"chunk_0_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
            SeqIO.write(record, tmp_fasta, "fasta")
            tmp_fasta_path = tmp_fasta.name
            tmp_fastas.append(tmp_fasta_path)
        
        # If it is smaller than the chunk size, just add it
        chunk_jobs.append((record, 0, record.id, species, tmp_fasta_path))
    
    return chunk_jobs, tmp_fastas

def prepare_jobs_parallel(c: Config, n_processes=None):
    if n_processes is None:
        n_processes = min(c.max_threads, mp.cpu_count())
    
    # Read all records first
    records = list(SeqIO.parse(c.input_fasta, "fasta"))
    print(f"Processing {len(records)} sequence(s) using {n_processes} threads...")
    
    # Create a partial function with the config
    process_func = partial(process_single_record, config=c)
    
    # Combine all results
    all_chunk_jobs = []
    all_tmp_fastas = []
    
    # Process records in parallel with progress bar
    with mp.Pool(processes=n_processes) as pool:
        with tqdm(total=len(records), desc="Processing sequence(s)") as pbar:
            for result in pool.imap_unordered(process_func, records):
                chunk_jobs, tmp_fastas = result
                all_chunk_jobs.extend(chunk_jobs)
                all_tmp_fastas.extend(tmp_fastas)
                pbar.update()
    
    return all_chunk_jobs, all_tmp_fastas 
    
def parse_args():
    parser = argparse.ArgumentParser(
        description="Take a FASTA file and it through Blatdiver."
    )

    parser.add_argument(
        "-q", "--query", required = True, help="Path to the query FASTA file")
    parser.add_argument(
        "-o", "--output", required = True, help="Name of your output")
    parser.add_argument(
        "-d", "--database", required = True, help="Path to the blat database to use")
    parser.add_argument(
        "-tr", "--tree", required = True, help="The Time Tree of Life tree TimeTree_v5_Final.nwk")
    parser.add_argument(
        "-i", "--index", required = True, help="Index of sequences in your blat db and their species and locations")
    parser.add_argument(
        "-k", "--kraken", required = True, help="Path to Kraken2 database")
    parser.add_argument(
        "-t", "--threads", type=int, help="Specify the number of threads you want to use. Highly recommended to use multiple threads to decrease time. Default: 1")
    parser.add_argument(
        "-c", "--chunk", type=int, help="If you want to specify the size of the how we will split the sequence to feed it to blat, Default: 100000")
    parser.add_argument(
        "-r", "--remove", type = int, help = "If you don't want the program to clean up its files, set to 0. Default: 1 (True)")
    parser.add_argument(
        "-s", "--species", type = str, help = "If you know the species of your sequence, or all sequences, you can define it here (make sure spaces are replaced with '_'. If you are feeding in multiple species within one query, this will assume they are all the same species")
    parser.add_argument(
        "-minScore", type = int, help = "Is the minimum alignment score threshold, sets blat's parameter -minScore, default 50")
    parser.add_argument(
        "-minIdentity", type = int, help = "Is the threshold of the minimum percent identity a hit must have. Default: 0. Percent identity = ( match / Q_end - Q_start )*100")
    parser.add_argument(
        "-overlap_filter", type = int, help = "If you want to do additional filtering with the overlap filter, set to 1. Default: 0 (Filter explained in https://github.com/graceoualline/blatdiver/)")
    parser.add_argument(
        "-overlap_div_filter", type = int, help = "If you want to do additional filtering with the overlap and divergence filter, set to 1. Default: 0 (Filter explained in https://github.com/graceoualline/blatdiver/)")

    return parser.parse_args()

def main():
    args = parse_args()
    c = Config(
    input_fasta=args.query,
    chunk_size=args.chunk if args.chunk else 100000,
    output_dir=args.output,
    database = Path(args.database),
    index = args.index,
    max_threads=args.threads if args.threads else 1,
    kraken_db = args.kraken,
    tree = Tree(args.tree),
    remove = bool(args.remove) if args.remove == 0 else True,
    species = args.species if args.species else None,
    minScore = args.minScore if args.minScore else 30,
    minIdentity = args.minIdentity if args.minIdentity else 0,
    overlap_filter = bool(args.overlap_filter) if args.overlap_filter==1 else False,
    overlap_div_filter = bool(args.overlap_div_filter) if args.overlap_div_filter==1 else False
    )

    if not os.path.exists(c.output_dir):
        os.makedirs(c.output_dir)

    print("Start time:", datetime.now())

    #check and ensure the database is sound
    blat_files = [Path(f) for f in os.listdir(c.database) if f.endswith(".2bit")]
    blat_files.sort()
    ooc_files = [Path(f) for f in os.listdir(c.database) if f.endswith(".ooc")]
    ooc_files.sort()
    if len(blat_files) != len(ooc_files):
        raise ValueError(
            f"Mismatch in file counts:\n"
            f"- Found {len(blat_files)} BLAT file(s)\n"
            f"- Found {len(ooc_files)} OOC file(s)\n"
            "Please ensure every BLAT file has a matching OOC file."
        )
    
    # go through and build jobs
    chunk_jobs, tmp_fastas = prepare_jobs_parallel(c)
    
    print(f"Prepared {len(chunk_jobs)} jobs(s) to process.")

    #this will parallelize it
    #run the function run_blat_on_chunk on all different threads
    print(f"Running on {c.max_threads} threads")
    run_blat_on_chunk(chunk_jobs, blat_files, ooc_files, c)
    print("All BLAT jobs finished.")

    print("Now filtering Blat Output")
    run_all_filters(chunk_jobs, c)
    print("All filtering complete")

    print("Creating final result files")
    combine_all_results(c)

    #remove all of the tmp fastas that were created
    for tmp in tmp_fastas:
        os.remove(tmp)
    print("BLATDIVER FINISHED")

    print("Endtime time:", datetime.now())
    
main()