#main place to run the blat command
from filter_blat import *
import sys
from ete3 import Tree
import subprocess
import os
from pathlib import Path
import glob
from seq_id_index import *
from filter_blat_multi import *
from find_overlap_and_div import *


#have the number of db be 1000 for now, fix later for public use
# same for number of threads
#later do the split file thing so people without large CPU can still get fast blat
def run_blat(directory, query, output, blat_dir, num_blat_db = 1000, num_threads = 32):

    blat_dir = Path(blat_dir)  # Ensure blat_dir is a Path object
    #print("directory in run blat", directory)
    #print("output in run blat", output)
   
    for i in range(1, num_blat_db + 1):
        #output_path = os.path.join(directory, f"{output}_part_{i}.psl")
        output_path = f"{output}_part_{i}.psl"

        # Construct the command
        command = (
            f"blat {blat_dir / f'gtdb1k_part_{i}_output.2bit'} {query} "
            f"-ooc={blat_dir / f'part{i}_1k_11.ooc'} -tileSize=11 "
            f"{output_path} -q=dna -t=dna"
        )
        print("running", command)

        # Run the command
        subprocess.run(command, shell=True, cwd=directory)
        #subprocess.run(command, shell=True)

def get_q_species(input_file):
    #implement this later using kraken
    return "Enterococcus faecalis"

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    # if q_species = unk, then run the get query species
    if len(sys.argv) != 7:
        print("Usage: python3 blat_main.py <input query file.fasta> <output file name> <gtdb seq species index .pkl> <div tree> <blat_db_dir> <query species (joined with '_')>")
        sys.exit(1)

    # Usage example:
    input_file = sys.argv[1]
    output_name = sys.argv[2]
    #print("output_name in blat_main", output_name)
    kraken = sys.argv[3]
    tree = Tree(sys.argv[4])
    blat_dir = sys.argv[5]
    q_species = sys.argv[6]
    #TODO implement this
    if q_species == "unk": q_spec = get_q_species(input_file)

    #make directory where all blat output will go
    directory = output_name + "_blat_output"
    #print("directory", directory)
    os.makedirs(directory, exist_ok=True)  # Creates the directory if it doesn't exist

    # Step 2: Run all 1000 blat runs in that directory
    #remove the previous directories of the output name so it doesn't get confused
    out = output_name.split("/")[-1]
    run_blat(directory, input_file, out, blat_dir)

    #then filter all blat results and add divergence or ani
    input_files = glob.glob(os.path.join(directory, "*.psl"))
    output_file = output_name + "_blatdiver_output.tsv"

    kraken = load_hash_table(kraken)
    #do initial filter to keep everyone that is divergent
    filter_blat_multi(input_files, output_file, q_species, kraken, tree, input_file, blat_dir)

    #do second filter for overlap of chunks that map to divergently different species
    output_filtered_file = output_name + "_blatdiver_output_overlap_div_filtered.tsv"
    rows = compress(output_file)
    find_overlap_and_div(rows, output_filtered_file, tree, blat_dir, kraken)
 
#python3 blat_pipeline/blat_main.py /usr1/gouallin/blastdiver/test_plasmids/one_lucky_plasmid.fasta test /usr1/shared/all_gtdb_id_and_kraken_species.txt /usr1/shared/TimeTree_v5_Final.nwk /usr1/shared/gtdb_split_2bit_1k/

#testing to make sure ani function works

#python3 blat_pipeline/blat_main.py /usr1/gouallin/blat/rand_hidden_plasmids_150kbp/NC_016607.1,2488,25624,CP073938.1.fasta NC_016607.1,2488,25624,CP073938.1 /usr1/shared/all_gtdb_id_and_kraken_species_sorted.txt /usr1/shared/TimeTree_v5_Final.nwk /usr1/shared/gtdb_split_2bit_1k/ Escherichia_coli
'''
input_files = glob.glob(os.path.join("NZ_CP013043.1,12719,71166,CP137225.1_blat_output/", "*.psl"))
output_filtered_file = "blatdiver_out_NZ_CP013043.1_troubleshoot.tsv"
q_species = "Proteus_mirabilis"
kraken = load_hash_table("/usr1/shared/all_gtdb_seq_species_index.pkl")
tree = Tree("/usr1/shared/TimeTree_v5_Final.nwk")
q_seq = "/usr1/gouallin/blat/rand_hidden_plasmids_150kbp/NZ_CP013043.1,12719,71166,CP137225.1.fasta"
blat_dir =  "/usr1/shared/gtdb_split_2bit_1k/"
filter_blat_multi(input_files, output_filtered_file, q_species, kraken, tree, q_seq, blat_dir)
'''
