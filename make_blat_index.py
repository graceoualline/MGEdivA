# VERY SLOW DO NOT USE
# I will intake the kraken file with all seq ids in the first column and all species in the second
# then I will find the file name of the 2bit file where the seq id is located
# and add it to the 3rd column

import os 
import subprocess

def process_value(ref_id):
    blat_db = "/usr1/shared/gtdb_split_2bit_1k/"
    
    # Step 1: Identify which .2bit file contains the reference sequence
    ref_2bit_file = None
    #Eventuallt replace this search with just an index of ref id and where its located
    print(f"searching for {ref_id}")
    for file in os.listdir(blat_db):
        if file.endswith(".2bit"):
            #print("searching file", file)
            file_path = os.path.join(blat_db, file)
            result = subprocess.run(["twoBitInfo", file_path, "stdout"], capture_output=True, text=True)
            if ref_id in result.stdout:
                ref_2bit_file = file_path
                print(f"{ref_id} found in", ref_2bit_file)
                break
    if not ref_2bit_file:
        print(f"Error: Reference ID {ref_id} not found in {blat_db}")
        return None
    ref_file = ref_2bit_file.split("/")
    return ref_file[-1]

def main():
    kraken_file = "/usr1/shared/all_gtdb_id_and_kraken_species.txt"
    output_file = "all_gtdb_id_species_location.tsv"  # Replace with your desired output file


    with open(kraken_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            columns = line.strip().split("\t")  # Adjust delimiter if needed
            seq_id = columns[0]
            seq_file = process_value(seq_id)  # Run first column through function
            out_line = "\t".join(columns) + f"\t{seq_file}\n"
            outfile.write(out_line)  # Store line in the list


    print(f"Output written to {output_file}")
