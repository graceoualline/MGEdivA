# blatdiver for optimized threading
# will run all blat db on different threads
# then run filters on all results on different threads
# combine all results, filters and blat output into one file
# remove individual chunk files (if sequence over 100kbp)
# or if there are differerent sequences


import os
import tempfile
import subprocess
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from blat_main import *
from datetime import datetime
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
    minScore: int # Is the minimum alignment score threshold, sets blat's parameter -minScore, default 50
    minIdentity: int # Is the threshold of the minimum percent identity a hit must have. Default: 95. Percent identity = ( match / Q_end - Q_start )*100

def combine_and_cleanup_psl_files(args):
    header_lines = 5
    output_chunk_dir, combined_path, c = args

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
                if c.remove: os.remove(fname)
        if c.remove: os.rmdir(output_chunk_dir)

        #print(f"Combined {len(input_files)} files into {combined_path}.")

def run_div_filter(job):
    chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = job
    try:
        #run the first round of filtering on the raw blat data
        filter_blat(all_blat_raw, output_file, species, c.index, c.tree, tmp_fasta_path, c.database)
    except subprocess.CalledProcessError as e:
        print(f"Error running filter_blat on chunk {output_file}:\n{e.stderr.decode().strip()}")
    return 1

def run_specified_filter(job):
    """Dispatch the correct function based on job type."""
    job_type, args = job
    try:
        if job_type == "no_gaps":
            to_filter_file, gap_file = args
            remove_large_gaps(to_filter_file, gap_file)
        elif job_type == "overlap":
            infile, outfile, c = args
            rows = compress(infile)
            find_overlap(rows, outfile, c.database, c.index)
        elif job_type == "overlap_div":
            to_filter_file, overlap_div_file, c = args
            rows = compress(to_filter_file)
            find_overlap_and_div(rows, overlap_div_file, c.tree, c.database, c.index)
    except subprocess.CalledProcessError as e:
        print(f"[{job_type}] Error on {args[1] if len(args) > 1 else 'unknown'}:\n{e.stderr.decode().strip()}")
    return 1

def run_no_gaps(job):
    to_filter_file, gap_file = job
    try:
        #run the first round of filtering on the blatdiver output
        remove_large_gaps(to_filter_file, gap_file)
    except subprocess.CalledProcessError as e:
        print(f"Error running no_gaps on {gap_file}:\n{e.stderr.decode().strip()}")
    return 1

def run_overlap(job):
    to_filter_file, overlap_output_file, c = job
    try:
        #run the first round of filtering on the blatdiver output
        rows = compress(to_filter_file)
        find_overlap(rows, overlap_output_file, c.database, c.index)
    except subprocess.CalledProcessError as e:
        print(f"Error running overlap on {overlap_output_file}\n{e.stderr.decode().strip()}")
    return 1

def run_overlap_div(job):
    to_filter_file, overlap_div_file, c = job
    try:
        #run the first round of filtering on the blatdiver output
        rows = compress(to_filter_file)
        find_overlap_and_div(rows, overlap_div_file, c.tree, c.database, c.index)
    except subprocess.CalledProcessError as e:
        print(f"Error running overlap div on {overlap_div_file}:\n{e.stderr.decode().strip()}")
    return 1

def run_all_filters(chunk_jobs, c: Config):
    blat_combine_jobs = []
    jobs = []
    for args in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = args
        output_chunk_dir = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}")
        
        #make one big file that combines all of the blat outputs
        all_blat_raw = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if not(os.path.exists(all_blat_raw)): blat_combine_jobs.append((output_chunk_dir, all_blat_raw, c))
        
        #then filter all blat results and add divergence or ani
        output_file = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blatdiver_output.tsv")
        jobs.append((chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c))

    #first, combine all blat files
    if len(blat_combine_jobs) > 0:
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(blat_combine_jobs), desc="Combining blat data") as pbar:
                for _ in pool.imap_unordered(combine_and_cleanup_psl_files, blat_combine_jobs):
                    pbar.update()
    print("Finished combining blat data.")
    
    
    #run divergence filter
    div_jobs = []
    for j in jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = j
        if os.path.exists(output_file): print(f"Skipping {output_file}, already exists.")
        else: div_jobs.append(j)
    if len(div_jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(div_jobs), desc="Running chunks through first divergence filter") as pbar:
                for _ in pool.imap_unordered(run_div_filter, div_jobs):
                    pbar.update()
    print("Finished running divergence filter.")

    #do gap jobs
    gap_jobs = []
    for j in jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = j
        gaps = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_no_gaps.tsv")
        if os.path.exists(gaps): print(f"Skipping {gaps}, already exists.")
        else: gap_jobs.append(("no_gaps", (output_file, gaps)))
    if len(gap_jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(gap_jobs), desc="Running chunks through no gaps filter") as pbar:
                for _ in pool.imap_unordered(run_specified_filter, gap_jobs):
                    pbar.update()
    print("Finished running no gaps filter")


    #run overlap, overlap_div, overlap_gav, overlap_gap_div filters
    filter_jobs = []
    for j in jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, blatdiver_output_file, c = j

        #overlap filter
        overlap = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_overlap.tsv")
        if os.path.exists(overlap): print(f"Skipping {overlap}, already exists.")
        else: filter_jobs.append(("overlap", (blatdiver_output_file, overlap, c)))

        # overlap div filter
        overlap_div = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_filtered_overlap_div.tsv")
        if os.path.exists(overlap_div): print(f"Skipping {overlap_div}, already exists.")
        else: filter_jobs.append(("overlap_div", (blatdiver_output_file, overlap_div, c)))

        chunk_record, idx, original_id, species, tmp_fasta_path, all_blat_raw, output_file, c = j
        gap = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_no_gaps.tsv")

        #overlap_gap filter
        overlap_gap = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_overlap_gap.tsv")
        if os.path.exists(overlap_gap): print(f"Skipping {overlap_gap}, already exists.")
        else: filter_jobs.append(("overlap", (gap, overlap_gap, c)))
        # overlap div gap filter
        overlap_div_gap = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_no_gap_overlap_div.tsv")
        if os.path.exists(overlap_div_gap): print(f"Skipping {overlap_div_gap}, already exists.")
        else: filter_jobs.append(("overlap_div", (gap, overlap_div_gap, c)))
    
    #run this first round of filtering
    if len(filter_jobs) > 0:
        with Pool(processes=c.max_threads) as pool:
            #Progress bars
            with tqdm(total=len(filter_jobs), desc="Running chunks through overlap, overlap div, overlap gap, and overlap_gap_div filters") as pbar:
                for _ in pool.imap_unordered(run_specified_filter, filter_jobs):
                    pbar.update()
    print("All filters complete")
    
def run_blat(args):
    blat_dir, blat_file, ooc_file, query, output_path, c = args
    # Construct the command
    command = (
        f"blat {blat_dir / blat_file} {query} "
        f"-ooc={blat_dir / ooc_file} -tileSize=11 -minScore={c.minScore}"
        f"{output_path} -q=dna -t=dna"
    )
    try:
        #run the command
        if os.path.exists(output_path): print(f"{output_path} skipped, already exists")
        else: 
            #print("Running", command)
            subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAT on {output_path}:\n{e.stderr.decode().strip()}")
    return 1

def run_blat_on_chunk(chunk_jobs, blat_files, ooc_files, c: Config):
    num_blat_db = len(blat_files)

    jobs = []
    # for each chunk, make all of the blat search jobs
    for chunk in chunk_jobs:
        chunk_record, idx, original_id, species, tmp_fasta_path = chunk
    
        # make a temporary file that isolates the fasta on its ow
        #get the name and path for the output file
        all_blat_output = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}_blat_output.tsv")
        if os.path.exists(all_blat_output): print(f"Skipping BLAT on chunk_{idx}_{original_id}, chunk_{idx}_{original_id}_blat_output.tsv already exists.")
        else:
            output_chunk = os.path.join(c.output_dir, f"chunk_{idx}_{original_id}")
            if not os.path.exists(output_chunk): os.makedirs(output_chunk)
            
            for i in range(num_blat_db):
                output_path = f"chunk_{idx}_{original_id}_part_{i}.psl"
                output_path = os.path.join(output_chunk, output_path)
                jobs.append((c.database, blat_files[i], ooc_files[i], tmp_fasta_path, output_path, c))
    if len(jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(jobs), desc="Running Sequences through BLAT") as pbar:
                for _ in pool.imap_unordered(run_blat, jobs):
                    pbar.update()

def run_combine_results(job):
    output_dir, chunk_size, to_combine, output_path, c = job
    try: 
        to_remove = adjust_and_merge_tsvs(output_dir, chunk_size, output_path, to_combine)
        if c.remove: 
            for f in to_remove: os.remove(f)
    except subprocess.CalledProcessError as e:
        print(f"Error combining results {output_path}: {e}")

def combine_all_results(c: Config):
    #once all jobs are finished, combine all of the filters into one place:
    parent_dir = os.path.dirname(c.output_dir)
    output_name = os.path.basename(c.output_dir)
    #python3 combine_blatdiver_output.py <blatdiver_output_directory> <chunk_size> <which files to combine ex. *blatdiver_output.tsv> <output_file>
    input_output = [("blat_output.tsv", os.path.join(c.output_dir, f"{output_name}_blat_results.tsv")),
                    ("filtered_overlap_div.tsv", os.path.join(c.output_dir, f"{output_name}_overlap_div.tsv")), 
                    ("blatdiver_output.tsv", os.path.join(c.output_dir, f"{output_name}_blatdiver_output.tsv")),
                    ("no_gap_overlap_div.tsv", os.path.join(c.output_dir, f"{output_name}_no_gaps_overlap_div.tsv")),
                    ("no_gaps.tsv", os.path.join(c.output_dir, f"{output_name}_no_gaps.tsv")), 
                    ("overlap_gap.tsv", os.path.join(c.output_dir, f"{output_name}_overlap_gap.tsv")),
                    ("overlap.tsv", os.path.join(c.output_dir, f"{output_name}_overlap.tsv"))]

    jobs = []
    for to_combine, output in input_output:
        output_path = os.path.join(parent_dir, output)
        if not os.path.exists(output_path): jobs.append((c.output_dir, c.chunk_size, to_combine, output_path, c))
    if len(jobs) > 0:
        #print whatever command you are running
        with Pool(processes=c.max_threads) as pool:
            #Progress bar
            with tqdm(total=len(jobs), desc="Combining all results") as pbar:
                for _ in pool.imap_unordered(run_combine_results, jobs):
                    pbar.update()

def prepare_jobs(c:Config):
    chunk_jobs = []
    tmp_fastas = []
    for record in SeqIO.parse(c.input_fasta, "fasta"):
        #get the length of the sequence
        seq_len = len(record.seq)
        #get the species of the sequence
        if c.species == None:
            try:
                with tempfile.NamedTemporaryFile(prefix=f"{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                    SeqIO.write(record, tmp_fasta, "fasta")
                    tmp_fasta_path = tmp_fasta.name
                species = get_q_species(tmp_fasta_path, c.kraken_db)
                os.remove(tmp_fasta_path)
            except:
                print("Issue extracting species for", record[:0])
                species = "unclassified"
        else: species = c.species
        print("Species:", species)
        #if the sequence length is greater than the chunk size we want
        if seq_len > c.chunk_size:
            for i in range(0, seq_len, c.chunk_size):
                #then we chunk it by size
                chunk_seq = record.seq[i:i+c.chunk_size]
                chunk_record = record[:0]  # copy header
                #set the sequence to be this chunk
                chunk_record.seq = chunk_seq
                #set the id to be the chunk number
                idx = i // c.chunk_size
                chunk_record.id = f"{record.id}_chunk_{idx}"
                chunk_record.description = f"{record.description} chunk {idx}"
                #get a tmp file of the sequence
                with tempfile.NamedTemporaryFile(prefix=f"chunk_{idx}_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                    SeqIO.write(chunk_record, tmp_fasta, "fasta")
                    tmp_fasta_path = tmp_fasta.name
                    tmp_fastas.append(tmp_fasta_path)
                #record the actual sequence and its record, its chunk number, and its id
                chunk_jobs.append((chunk_record, idx, record.id, species, tmp_fasta_path))
        else:
            #get a tmp file of the sequence
            with tempfile.NamedTemporaryFile(prefix=f"chunk_0_{record.id}_tmp", mode="w", delete=False, suffix=".fa") as tmp_fasta:
                SeqIO.write(record, tmp_fasta, "fasta")
                tmp_fasta_path = tmp_fasta.name
                tmp_fastas.append(tmp_fasta_path)
            #if it is smaller than the chunk size, just add it
            chunk_jobs.append((record, 0, record.id, species, tmp_fasta_path)) 
    return chunk_jobs, tmp_fastas   
    
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
        "-minIdentity", type = int, help = "Is the threshold of the minimum percent identity a hit must have. Default: 95. Percent identity = ( match / Q_end - Q_start )*100")

    return parser.parse_args()

def main():
    args = parse_args()
    c = Config(
    input_fasta=args.query,
    chunk_size=args.chunk if args.chunk else 100000,
    output_dir=args.output,
    database=Path(args.database),
    index=load_hash_table(args.index),
    max_threads=args.threads if args.threads else 1,
    kraken_db=args.kraken,
    tree= Tree(args.tree),
    remove = bool(args.remove) if args.remove == 0 else True,
    species = args.species if args.species else None,
    minScore = args.minScore if args.minScore else 50,
    minIdentity = args.minIdentity if args.minIdentity else 95
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
    chunk_jobs, tmp_fastas = prepare_jobs(c)
    
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