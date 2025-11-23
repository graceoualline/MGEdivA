#this will intake a psl blat file and only output the results where there is divergence between the two species
# will also add a column that will show the divergence time between the query and reference species, and those two species

#assuming we figured out the query species in the beginning using Kraken

import pandas as pd
import sys
from ete3 import Tree
import subprocess
#import glob
import os
import tempfile
#import bisect
from build_database_index import *
from Bio import SeqIO
#from collections import defaultdict
#from twobitreader import TwoBitFile

# def extract_2bit_file_batch(file_path_2bit, seq_ids):
#     tb = TwoBitFile(file_path_2bit)
#     seq_file_dict = dict()

#     for name in seq_ids:
#         if name in tb:
#             with tempfile.NamedTemporaryFile(prefix=f"refs{name}", suffix=".fa", delete=False, mode='w') as temp_fa:
#                     temp_path = temp_fa.name
#             seq = tb[name]
#             seq_file_dict[name] = temp_path
#             with open(temp_path, "w") as out:
#                 out.write(f">{name}\n{seq}\n")

#     return seq_file_dict

def batch_extract_sequences(ref_ids, index, blat_db):
    """
    Extract multiple sequences at once, grouped by 2bit file.
    Returns dict of ref_id -> fasta_path
    """
    # Group ref_ids by their 2bit file
    # Extract all sequences from each 2bit file
    ref_to_fasta = {}
    fasta_to_ref = {}
    # for each 2bit file, extract these files
    i = 0
    for ref in ref_ids:
        i += 1
        print("file", i)
        try:
            #subprocess.run(cmd, check=True, capture_output=True)
            ref_to_fasta[ref] = extract_2bit_fasta(ref, index, blat_db)
            fasta_to_ref[ref_to_fasta[ref]] = ref
            
        except subprocess.CalledProcessError as e:
            print(f"Error extracting sequence {ref}: {e}")
            continue

    print("ref to fasta", len(ref_to_fasta))
    return ref_to_fasta, fasta_to_ref

def batch_calculate_ani(query_path, ref_paths_dict, paths_ref_dict):
    """
    Calculate ANI for all query-ref pairs using skani in batch mode.
    Returns dict of ref_id -> ANI value
    """
    if not ref_paths_dict:
        return {}
    # Create a temporary file listing all reference sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as ref_list:
        ref_list_path = ref_list.name
        for ref_path in ref_paths_dict.values():
            ref_list.write(f"{ref_path}\n")
    
    ani_results = {}
    
    try:
        # Run skani with query against all references
        print(f"skani dist {query_path} --rl {ref_list_path}")
        result = subprocess.run(
            ["skani", "dist", query_path, "--rl", ref_list_path],
            capture_output=True,
            text=True,
            check=True
        )
        
        # Parse output
        lines = result.stdout.strip().split("\n")
        with open("debugging_batch_skani.txt", "w") as f:
            for line in lines[1:]:  # Skip header
                f.write(f"{line}\n")
                parts = line.split()
                if len(parts) >= 3:
                    ref_path = parts[0]
                    ani_value = float(parts[2])
                    if ref_path in paths_ref_dict: ani_results[paths_ref_dict[ref_path]] = ani_value
        # Return -1 for all refs that failed
        for ref_id in ref_paths_dict.keys():
            if ref_id not in ani_results:
                ani_results[ref_id] = 0
                
    except (subprocess.CalledProcessError, ValueError) as e:
        print(f"Error in batch ANI calculation: {e}")
        with open("debugging_batch_skani_err.txt", "w") as f: f.write(e)
        input("Press enter to continue")
        # Return -1 for all refs that failed
        for ref_id in ref_paths_dict.keys():
            if ref_id not in ani_results:
                ani_results[ref_id] = "unk:skani_err"
    finally:
        os.remove(ref_list_path)
    
    return ani_results

# def get_path_from_name(leaf_name, tree):
#     if leaf_name == "NA": return "NA"
#     try: return tree & leaf_name
#     except: 
#         try: 
#             path = tree.get_path(leaf_name)
#             return path
#         except:
#             print("Unable to find path of", leaf_name)
#             return "NA"

#this is taken from PAReTT and altered
# def get_div(path1, path2, tree):
#     if path1 == "NA" or path2 == "NA":
#         return "unk:unable_to_find_ref_species_in_tree"
#     try:
#         distance = path1.get_distance(path2)
#         return distance / 2
#     except:
#         print("Unable to find the path betweens", path1, path2)
#         return "unk:unable_connect_two_paths_in_tree"


def calculate_distance(seq1, seq2):
    """
    Calculates the ANI distance between two sequences using skani.
    
    Parameters:
    - seq1 (str): Path to the query sequence (FASTA format).
    - seq2 (str): Path to the reference sequence (FASTA format).

    Returns:
    - float: ANI value (or 0 if not found).
    """
    #skani wont work if sequence under 500 bp
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
            #that means nothing was found, so we don't want this value
            return -1

    except (IndexError, ValueError, subprocess.CalledProcessError) as e:
        #print(f"Error running skani: {e}, if sequence was under 500 bp, skani will not work")
        return "unk:skani_err"  # Return 0 in case of errors

def extract_2bit_fasta(ref_id, index, blat_db):
    # Step 1: Identify which .2bit file contains the reference sequence
    ref_2bit_file = lookup_location(index, ref_id)
    if not ref_2bit_file:
        print(f"Error: Reference ID {ref_id} not found in {blat_db}")
        return None
    
    ref_2bit_file = os.path.join(blat_db, ref_2bit_file.lstrip("/"))
    
    # Step 2: Extract the reference sequence from the identified .2bit file
    with tempfile.NamedTemporaryFile(prefix = f"reference_{ref_id}", suffix=".fa", delete=False) as ref_fasta:
        ref_fasta_name = ref_fasta.name
        subprocess.run(["twoBitToFa", ref_2bit_file, "-seq=" + ref_id, ref_fasta_name], check=True)
    return ref_fasta_name

def find_ani(q_seq, ref_id, blat_db, index):
    """
    Finds the ANI between a query sequence and a reference sequence stored in a .2bit file.

    Parameters:
    - q_seq (str): Path to the query sequence (FASTA format).
    - ref_id (str): The reference sequence ID to extract.
    - blat_db (str): Path to the folder containing .2bit files.

    Returns:
    - float: ANI value (or None if reference is not found).
    """
    ref_fasta_name = extract_2bit_fasta(ref_id, index, blat_db)

    distance = calculate_distance(q_seq, ref_fasta_name)
    
    os.remove(ref_fasta_name)
    return distance

def check_div_cache(sp1, sp2, cache):
    if (sp1, sp2) in cache: return cache[(sp1, sp2)]
    elif (sp2, sp1) in cache: return cache[(sp2, sp1)]
    return None

def filter_blat(inf, outf, q_species, index, tree, q_seq, blat_db, minIdentity):
    div_cache = dict() # will be (species1, species2): divergence time
    #path_cache = dict() # will be species: name of it in the tree
    ani_cache = dict() # will be ref_id: ani
    #get a df or array that converts all of the columns fomr input file .psl into something readable

    # get the name of the q species that is in the tree
    q_og_species = q_species
    if q_species != "unclassified": q_species = tree.get_path_identifier(tree.get_path(q_species))
    
    lines_to_process = []
    ref_ids_needing_ani = set()

    with open(inf, 'r') as infile, open(outf, 'w') as outfile:
        
        outfile.write("match\tmismatch\trep. match\tN's\tQ gap count \tQ gap bases\tT gap count\tT gap bases\tstrand\tQ name\tQ size\tQ start\t Q end\tT name\tT size\tT start\tT end\tblock count\tblockSizes\tqStarts\ttStarts\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\n")
        # Process each line in the BLAT file
        for line in infile:
            columns = line.strip().split('\t')
            #skip the header lines
            if len(columns) != 21 or columns[0] == 'match':
                continue
            query_id = columns[9]  # Assuming query sequence ID is in the 10th column
            ref_id = columns[13]  # Assuming reference sequence ID is in the 14th column
            match = int(columns[0])
            Qstart = int(columns[11])
            Qend = int(columns[12])
            perIdent = round((match / (Qend-Qstart))*100, 5)

            #skip this line if it doesnt meet the threshold
            if perIdent < minIdentity:
                continue
            # Get the species of the reference sequence using Kraken and grep
            ref_species = lookup_species(index, ref_id)
    
            if ref_species is None:
                print("ERROR, your index mapping database genomes to their species and locations was not made correctly.")
                print(f"Unable to find the species for {ref_id} in {index}. Skipping.")
                continue  # Skip this line if no species found for reference sequence
            if q_species in ["unclassified", "NA"] or ref_species in ["unclassified", "NA"]: # or path1 == None:
                div = "unk:unclassified_species"
            else:
                div = check_div_cache(q_species, ref_species, div_cache)
                if div == None:
                    ref_leaf_name = lookup_tree_leaf_name(index, ref_id)
                    # Get the divergence time between query species and reference species
                    div = tree.divergence(q_species, ref_leaf_name)
                    div_cache[(q_species, ref_species)] = div
            if div == None:
                div = "unk:unable_to_find_ref_species_in_tree"

            ani = "NA" #ani is NA for now unless defined later
            if isinstance(div, (int, float)) and div < 1:
                continue  # Skip lines with less than 1 MYA
            elif isinstance(div, str): #if div is unk, then find ani
                lines_to_process.append({
                'columns': columns,
                'ref_id': ref_id,
                'perIdent': perIdent,
                'ref_species': ref_species,
                'div': div
                })
                ref_ids_needing_ani.add(ref_id)
                continue
                #if ref_id in ani_cache: ani = ani_cache[ref_id]
                #else: 
                #    ani = find_ani(q_seq, ref_id, blat_db, index) 
                #    ani_cache[ref_id] = ani
                #if isinstance(ani, (int, float)):
                #    if ani < 0 or ani >= 95:
                #        continue #skip lines with 95 or more ani since that means they are same species
                        # also skip lines where ani is unknown, since there is no way to know its relation to the query

            # write everyone that made it through the filters
            new_line = "\t".join(columns)
            outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{ref_species}\t{div}\t{ani}\n")
        
        # now go through everyone that needs ani taken
        if ref_ids_needing_ani:
            print(f"Extracting {len(ref_ids_needing_ani)} reference sequences...")
            ref_fastas, fastas_ref = batch_extract_sequences(ref_ids_needing_ani, index, blat_db)
            
            print(f"Calculating ANI for {len(ref_fastas)} sequences...")
            ani_results = batch_calculate_ani(q_seq, ref_fastas, fastas_ref)
            ani_cache.update(ani_results)
            
            # Cleanup temporary files
            for fasta_path in ref_fastas.values():
                try:
                    os.remove(fasta_path)
                except:
                    pass

            # then write them to file if they pass
            for line_data in lines_to_process:
                columns = line_data['columns']
                ref_id = line_data['ref_id']
                perIdent = line_data['perIdent']
                div = line_data['div']
                
                ani = "NA"
                ani = ani_cache.get(ref_id, "unk:skani_err")
                if isinstance(ani, (int, float)):
                    if ani < 0 or ani >= 95:
                        continue
                # Write the line with the added divergence time, query species, and reference species
                new_line = "\t".join(columns)
                outfile.write(new_line.strip() + f"\t{perIdent}\t{q_og_species}\t{lookup_species(index, ref_id)}\t{div}\t{ani}\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 9:
        print("Usage: python3 filter_blat.py <input psl file> <output .tsv file> <query species (make sure _ instead of space)> <gtdb seq species index .pkl> <div tree> <input sequence> <blat db> <minIdentity>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    q_species = sys.argv[3]
    index = sys.argv[4]
    tree = Divergence_Tree_Preprocessed(sys.argv[5])
    q_seq = sys.argv[6]
    blat_db = sys.argv[7]
    index = load_hash_table(index)
    minIdentity = float(sys.argv[8])
    filter_blat(input_file, output_file, q_species, index, tree, q_seq, blat_db, minIdentity)