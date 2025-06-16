#this will intake a psl blat file and only output the results where there is divergence between the two species
# will also add a column that will show the divergence time between the query and reference species, and those two species

#assuming we figured out the query species in the beginning using Kraken

import pandas as pd
import sys
#from phylorank import Tree
from ete3 import Tree
import subprocess
import glob
import os
import tempfile
import bisect
from seq_id_index import *

def get_path(sp, tree):
    sp = "_".join(sp.split(" "))
    sp = "'"+sp+"'"

    try:
        path = tree & sp
    except:
        first = sp.split("_")
        try:
            new = first[0] + "_" + first[1] + "'"
            #print("new1", new)
            path = tree & new
        except:
            try: 
                new = first[0] + "'"
                #print("new2", new)
                path = tree & new
            except:
                #if all else fails, find the first instance of the name within another name
                for node in tree.traverse("preorder"):
                    if first[0] in node.name:
                        return node
                print("sp not found", sp)
                return None
    return path

#this is taken from PAReTT and altered
def get_div(path1, path2, tree):
    #print(sp1, sp2)
    #path1 = get_path(sp1, tree)
    #path2 = get_path(sp2, tree)
    #print(sp1, sp2)
    if path1 == None or path2 == None:
        return None

    distance = path1.get_distance(path2)
    
    return distance

def process_kraken(hash_table, seq_id):
    """
    Uses binary search to find the species of a given sequence ID from the Kraken output.
    Returns the species name, or None if not found.
    """
    # Use binary search to find the sequence ID in the Kraken file
    #print("seq_id", seq_id)
    try:
        species = lookup_species(hash_table, seq_id)
        return species
    except Exception as e:
        print(f"Error when finding {seq_id}: {e}")
        return "unclassified"

def calculate_distance(seq1, seq2):
    """
    Calculates the ANI distance between two sequences using skani.
    
    Parameters:
    - seq1 (str): Path to the query sequence (FASTA format).
    - seq2 (str): Path to the reference sequence (FASTA format).

    Returns:
    - float: ANI value (or 0 if not found).
    """
    try:
        # Run the skani command to calculate distance
        result = subprocess.run(
            ["skani", "dist", seq1, seq2],
            capture_output=True,
            text=True,
            check=True
        )
        # Capture the output and process it to extract the ANI value
        output = result.stdout
        lines = output.strip().split("\n")
        
        # Assume the second line contains the data, split by whitespace to extract the ANI
        if len(lines) > 1:
            ani_value = lines[1].split()[2]
            return float(ani_value)
        else:
            #raise ValueError("Unexpected output format from skani.")
            #that means nothing was found, so ANI is 0
            return 0
    except (IndexError, ValueError, subprocess.CalledProcessError) as e:
        print(f"Error running skani: {e}")
        return 0  # Return 0 in case of errors


def extract_2bit_fasta(ref_id, kraken, blat_db):
    # Step 1: Identify which .2bit file contains the reference sequence
    ref_2bit_file = lookup_location(kraken, ref_id)
    if not ref_2bit_file:
        print(f"Error: Reference ID {ref_id} not found in {blat_db}")
        return None
    
    ref_2bit_file = blat_db + "/" + ref_2bit_file
    
    # Step 2: Extract the reference sequence from the identified .2bit file
    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as ref_fasta:
        ref_fasta_name = ref_fasta.name
        subprocess.run(["twoBitToFa", ref_2bit_file, "-seq=" + ref_id, ref_fasta_name], check=True)
        #print(ref_fasta_name)
        #input("Press Enter after checking the file...") #pause the program so I can check the contents of the ref_file
    return ref_fasta_name

def find_ani(q_seq, ref_id, blat_db, kraken):
    """
    Finds the ANI between a query sequence and a reference sequence stored in a .2bit file.

    Parameters:
    - q_seq (str): Path to the query sequence (FASTA format).
    - ref_id (str): The reference sequence ID to extract.
    - blat_db (str): Path to the folder containing .2bit files.

    Returns:
    - float: ANI value (or None if reference is not found).
    """
    ref_fasta_name = extract_2bit_fasta(ref_id, kraken, blat_db)
        
    distance = calculate_distance(q_seq, ref_fasta_name)
    
    os.remove(ref_fasta_name)
    return distance

def filter_blat(inf, outf, q_species, kraken, tree, q_seq, blat_db):
    #get a df or array that converts all of the columns fomr input file .psl into something readable
    if q_species != "unclassified": path1 = get_path(q_species, tree)
    #write header of input file to putput file
    with open(inf, 'r') as infile, open(outf, 'w') as outfile:
        # Write a new header to file that is tabulated and has our new guys on the end
        #we only care about query, sequence
        #outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \t Q gap bases\tT gap count \tT gap bases\tstrand\tQ name\tQ size\tQ start\t Q end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tQuery_Species\tReference_Species\tDivergence_Time\n")
        outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \t Q gap bases\tT gap count \tT gap bases\tstrand\tQ name\tQ size\tQ start\t Q end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tQuery_Species\tReference_Species\tDivergence_Time\tANI_of_whole_seqs(only_if_div=unk)\n")
        species_path_cache = {}
        # Process each line in the BLAT file
        for line in infile:
            columns = line.strip().split('\t')
            #print("columns", columns)
            #skip the header lines
            if len(columns) != 21 or columns[0] == 'match':
                continue
            query_id = columns[9]  # Assuming query sequence ID is in the 10th column
            ref_id = columns[13]  # Assuming reference sequence ID is in the 14th column
            
            # Get the species of the reference sequence using Kraken and grep
            ref_species = process_kraken(kraken, ref_id)
    
            
            if ref_species is None:
                print("ERROR, impossible that we cannot find a sequence that is in the GTDB")
                assert(False)
                continue  # Skip this line if no species found for reference sequence
            if q_species == "unclassified" or ref_species == "unclassified" or path1 == None:
                div = "unk:unclassified_species"
            else:
                if ref_species not in species_path_cache:
                    species_path_cache[ref_species] = get_path(ref_species, tree)
            # Get the divergence time between query species and reference species
                div = get_div(path1, species_path_cache[ref_species], tree)
            print("div", div, q_species, ref_species)
            if div == None:
                div = "unk:unable_to_find_ref_species_in_tree"
                
            if div == 0:
                continue  # Skip lines with zero divergence time
            elif "unk" in div:
                print("div unknown:", ani)
                ani = find_ani(q_seq, ref_id, blat_db, kraken) 
                print("ani calculated:", ani)
                if ani >= 95:
                    continue #skip lines with 95 or more ani since that means they are same species
            else:
                ani = "N/A"
            
            # Write the line with the added divergence time, query species, and reference species
            new_line = columns #columns[9:]
            new_line = "\t".join(new_line)
            outfile.write(new_line.strip() + f"\t{q_species}\t{ref_species}\t{div}\n")


if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 8:
        print("Usage: python3 filter_blat.py <input psl file> <output .tsv file> <query species (make sure _ instead of space)> <gtdb seq species index .pkl> <div tree> <input sequence> <blat db>")
        sys.exit(1)

    # Usage example:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    q_species = sys.argv[3]
    kraken = sys.argv[4]
    tree = Tree(sys.argv[5])
    q_seq = sys.argv[6]
    blat_db = sys.argv[7]
    kraken = load_hash_table(kraken)
    filter_blat(input_file, output_file, q_species, kraken, tree, q_seq, blat_db)
    

#kraken = /usr1/shared/all_gtdb_id_and_kraken_species.txt
#tree = /usr1/shared/TimeTree_v5_Final.nwk

#in pwd /usr1/gouallin/blastdiver/blat
#python3 filter_blat.py blat_test_output100.psl filtered_blat_100.psl Klebsiella_pneumoniae_subsp._pneumoniae /usr1/shared/all_gtdb_id_and_kraken_species.txt /usr1/shared/TimeTree_v5_Final.nwk