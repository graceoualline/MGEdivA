#this will intake a kraken output file
# and return a .txt list of all sequence id's and species that were found
import tempfile
import os
import subprocess
import sys

def process_all_gtdb_kraken(kraken_file):
    # Use grep to find the sequence ID in the Kraken file
    all_species = set()
    with open(kraken_file, "r") as infile:
        # Extract the species name from the grep output
        # The species name is in the 3rd column (assuming Kraken output format is consistent)

        for line in infile:
            id = line.split('\t')[1]
            species = line.split('\t')[2].split('(')[0].strip()
            #print(id, species)
            all_species.add((id, species))
    return all_species

def get_sp_from_kraken(kraken_output_file):
# Use grep to find the sequence ID in the Kraken file
    with open(kraken_output_file, "r") as infile:
        # Extract the species name from the grep output
        # The species name is in the 3rd column (assuming Kraken output format is consistent)
        lines = infile.readlines()
        if not lines:
            print("Kraken output is empty.")
            return "unclassified"

        last_line = lines[-1]
        last_line = last_line.split()
        species = last_line[5:]
        species = "_".join(species)
        #print("Last line:", last_line)
        print(species)
        
    return species


def get_q_species(fasta_file, kraken_db):
    with tempfile.NamedTemporaryFile(prefix="kraken_output", delete=False, mode='w+') as tmp:
        tmp_filename = tmp.name
    command = [
        "kraken2", "--db", kraken_db, "--threads", "1",
        "--report", tmp_filename, fasta_file
    ]
    print("Finding query species")
    print(" ".join(command))
    # Run the command
    subprocess.run(command, check = True)
    print("Done running kraken")
    species = get_sp_from_kraken(tmp_filename)

    if os.path.exists(tmp_filename):
        os.remove(tmp_filename)
    
    return species
    #except:
    #    print("Error finding species of the query")
    #    if os.path.exists(tmp_filename):
    #        os.remove(tmp_filename)

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    # if q_species = unk, then run the get query species
    if len(sys.argv) != 3:
        print("Usage: python3 extract_species_from_kraken.py <query fasta file> <kraken db> ")
        sys.exit(1)
    q_seq = sys.argv[1]   
    kraken_db = sys.argv[2]    
    get_q_species(q_seq, kraken_db)
#db = "/usr1/shared/kraken2_custom_db"
#file = "/usr1/gouallin/random_ctrl/full_genomes/full_rand_1/CP003077.2_random.fasta"
#get_q_species(file, db)
'''
kraken_file = "kraken_all_gtdb_output.txt"
all_species = process_kraken(kraken_file)
print("num of ids and species", len(all_species))
output_file = "all_gtdb_id_and_kraken_species.txt"

with open(output_file, 'w') as outfile:
    for id, species in all_species:
        outfile.write(f"{id}\t{species}\n")'''