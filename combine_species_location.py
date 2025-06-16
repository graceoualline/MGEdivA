# File paths (you can change these)
species_file = "/usr1/shared/all_gtdb_id_and_kraken_species.txt"         # Format: seq_id \t species
location_file = "sequence_id_to_2bitfile_mapping.txt"       # Format: seq_id \t file_location
output_file = "all_gtdb_seq_species_location.txt"    # Output file

# Step 1: Load species data into a dictionary
seq_to_species = {}
with open(species_file, "r") as f:
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            seq_id, species = parts[0], parts[1]
            seq_to_species[seq_id] = species

# Step 2: Open output and read location data, join with species
with open(output_file, "w") as fout, open(location_file, "r") as fin:
    for line in fin:
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            seq_id, file_location = parts[0], parts[1]
            print("finding", seq_id)
            species = seq_to_species.get(seq_id, "unknown_species")
            fout.write(f"{seq_id}\t{species}\t{file_location}\n")