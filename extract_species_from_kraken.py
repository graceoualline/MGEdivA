#this will intake a kraken output file
# and return a .txt list of all sequence id's and species that were found


def process_kraken(kraken_file):
    # Use grep to find the sequence ID in the Kraken file
    all_species = set()
    with open(kraken_file, "r") as infile:
        # Extract the species name from the grep output
        # The species name is in the 3rd column (assuming Kraken output format is consistent)

        for line in infile:
            id = line.split('\t')[1]
            species = line.split('\t')[2].split('(')[0].strip()
            print(id, species)
            all_species.add((id, species))
    return all_species

kraken_file = "kraken_all_gtdb_output.txt"
all_species = process_kraken(kraken_file)
print("num of ids and species", len(all_species))
output_file = "all_gtdb_id_and_kraken_species.txt"

with open(output_file, 'w') as outfile:
    for id, species in all_species:
        outfile.write(f"{id}\t{species}\n")