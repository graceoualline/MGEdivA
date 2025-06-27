#tihs will let me run filters on all of the outputs on different threads


import glob
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from multiprocessing import Pool

import sys
import os

def run_filter(args):
    # Adjust this command if needed
    filter_script, file_input, output_tail, other_info = args
    out = file_input.split("_")[:4]
    out = "_".join(out) + output_tail
    if os.path.exists(out):
        print(f"Skipping {out}, already exists.")
        return None  # Or you can return something like "Skipped"
    command = ["python3", filter_script, file_input, out] + other_info
    print("Running", command)
    return subprocess.run(command, capture_output=True, text=True)

def run_multi(file_tails, threads, filter, output_tail, other_info):
    files = glob.glob("*"+file_tails)  # Grab all matching files
    print(f"Found {len(files)} files.")

    jobs = [(filter, f, output_tail, other_info) for f in files]

    with Pool(processes=threads) as pool:
        pool.map(run_filter, jobs)

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python3 run_multi_files.py <file_tail ex *_blatdiver_output> <threads> <program or filter to use> <output_tail> <any other info needed for command>")
        sys.exit(1) 
    file_tail = sys.argv[1]
    threads = int(sys.argv[2])
    filter = sys.argv[3]
    output_tail = sys.argv[4]
    other_info = sys.argv[5:]
    #print(file_tail)
    #print(threads)
    #print(filter)
    #print(output_tail)
    #print(other_info

    run_multi(file_tail, threads, filter, output_tail, other_info)