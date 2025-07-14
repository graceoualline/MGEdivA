from filter_blat import *
import sys
from ete3 import Tree
import subprocess
import os
from pathlib import Path
import glob
from seq_id_index import *

def filter_blat_multi(infiles, outf, q_species, kraken, tree, q_seq, blat_db):
    """
    Processes multiple BLAT output files line by line, extracting relevant fields,
    performing species lookup and divergence calculation, and writing results
    to a single output file.
    
    Args:
        infiles (list): List of input BLAT .psl files.
        outf (str): Output file path.
        q_species (str): Query species name.
        kraken (str): Path to Kraken database.
        tree (str): Path to species divergence tree.
    """
    if q_species != "unclassified": path1 = get_path(q_species, tree)
    

    # Open output file once
    with open(outf, 'w') as outfile:
        # Write the header
        outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \t Q gap bases\tT gap count \tT gap bases\tstrand\tQ name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\n")
        species_path_cache = {}
        species_path_cache["unclassified"] = "NA"
        # Process each input file line by line
        for inf in infiles:
            #print("filtering", inf)
            with open(inf, 'r') as infile:
                for line in infile:
                    columns = line.strip().split('\t')

                    # Skip header lines
                    if len(columns) != 21 or columns[0] == 'match':
                        continue

                    query_id = columns[9]  # Query sequence ID
                    ref_id = columns[13]   # Reference sequence ID
                    
                    # Get the species of the reference sequence using Kraken and grep
                    #print("getting species of", ref_id)
                    #this is still a big bottle neck
                    ref_species = process_kraken(kraken, ref_id)
                    #print("ref species", ref_species)
                    
                    if ref_species is None:
                        print(f"ERROR: Reference sequence {ref_id} not found in GTDB.")
                        ref_species == "unclassified"
                        

                    # Determine divergence time
                    if q_species == "unclassified" or ref_species == "unclassified" or path1 == "NA":
                        div = "unk:unclassified_species"
                    else:
                        #print("getting div")
                        if ref_species not in species_path_cache:
                            #if species == "unclassified_Arthrobacter": print("filter_blat", outf)
                            species_path_cache[ref_species] = get_path(ref_species, tree)
                        # Get the divergence time between query species and reference species
                        div = get_div(path1, species_path_cache[ref_species], tree)
                        #print("div", div)

                    if div is None:
                        div = "unk:unable_to_find_ref_species_in_tree"

                    if type(div) == str and "unk" in div:
                        #print("div unknown:", div)
                        ani = find_ani(q_seq, ref_id, blat_db, kraken) #make find_ani later
                        #print("ani calculated:", ani)
                        if type(ani) != str and ani >= 95:
                            continue  # Skip lines with ANI >= 95 (same species)
                    # Skip those too closely related
                    elif div < 1:
                        continue
                    else:
                        ani = "N/A"
                    
                    # Construct and write the processed line
                    columns[-1] = columns[-1][:-2]
                    new_line = "\t".join(columns)#"\t".join(columns[9:])  # Extract relevant columns
                    new_line.replace("\n", "")
                    outfile.write(f"{new_line}\t{q_species}\t{ref_species}\t{div}\t{ani}\n")

            #print(f"Finished processing: {inf}")

    print(f"Processing complete. Results written to {outf}.")


if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    # if q_species = unk, then run the get query species
    if len(sys.argv) != 8:
        print("Usage: python3 filter_blat_multi.py <directory of blat output> <input query file.fasta> <output file name> <gtdb seq species index .pkl> <div tree> <blat_db_dir> <query species (joined with '_')>")
        sys.exit(1)

    # Usage example:
    directory = sys.argv[1]
    input_file = sys.argv[2]
    output_name = sys.argv[3]
    kraken = sys.argv[4]
    tree = Tree(sys.argv[5])
    blat_dir = sys.argv[6]
    q_species = sys.argv[7]
    #TODO implement this
    if q_species == "unk": q_spec = get_q_species(input_file)


    #then filter all blat results and add divergence or ani
    input_files = glob.glob(os.path.join(directory, "*.psl"))
    output_filtered_file = output_name + "_blatdiver_output.tsv"
    if os.path.exists(output_filtered_file):
        print(f"Skipping {output_filtered_file} , already exists.")
    else:
        kraken = load_hash_table(kraken)
        filter_blat_multi(input_files, output_filtered_file, q_species, kraken, tree, input_file, blat_dir)
 
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
