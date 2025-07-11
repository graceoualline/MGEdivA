# this should intake the folder that contains a bunch of split of 2bit files
# and it will return a .pkl file which is a hash table that essential works like
# index[sequence_id] = (sequence_species, sequence_location_in_the_database)
import os
import subprocess
import tempfile
import pickle

#first step, get the locations of all sequence ids
def get_seq_id_loc(two_bit_dir, input_files):
    seq_id_loc = set()
    
    for inputf in input_files:
        two_bit_path = f"{two_bit_dir}{inputf}"

        #this will make a temp file
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_path = temp_file.name

        # Run twoBitInfo to get sequence info
        try:
            print("running")
            print("twoBitInfo", two_bit_path, temp_path)
            result = subprocess.run(["twoBitInfo", two_bit_path, temp_path],
                                    capture_output = True,
                                    text = True,
                                    check = True)
        
            # Read the sequence IDs and write with file name
            with open(temp_path, "r") as f:
                for line in f:
                    seq_id = line.split()[0]
                    seq_id_loc.add((seq_id, inputf))
        finally:
            os.remove(temp_path)
            #assert(False)
    return seq_id_loc

def combine_seq_spec_loc(species_file, seq_id_loc):
    
    # Step 1: Load species data into a dictionary
    seq_to_species = {}
    with open(species_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                seq_id, species = parts[0], parts[1]
                seq_to_species[seq_id] = species

    # Step 2: Open output and read location data, join with species
    for seq_loc in seq_id_loc:
        #make sure there are two entries
        if len(seq_loc) >= 2:
            seq_id, file_location = seq_loc[0], seq_loc[1]
            species = seq_to_species.get(seq_id, "unknown_species")
            seq_to_species[seq_id] = (species, file_location)
    return seq_to_species

def main():
    two_bit_dir = "/usr1/shared/gtdb_2bil_split_2bit/"
    input_files = [f"split_{i}_output.2bit" for i in range(1, 138)]
    final_output = "gtdb_2bil_seq_id_species_loc_index.pkl"
    species_file = "/usr1/shared/all_gtdb_id_and_kraken_species.txt"

    #first step, get the locations of all sequence ids
    #used to be get_seq_ids.py
    seq_id_loc = get_seq_id_loc(two_bit_dir, input_files)
    # next, we use a file that connects seq_id and species determined by kraken
    # and add that to the index # use to be combine_species_location.py
    seq_species_loc = combine_seq_spec_loc(species_file, seq_id_loc)

    # # now make a loadable index from this new table of id, species, location
    with open(final_output, "wb") as f:
        pickle.dump(seq_species_loc, f)

    print(f"Hash index saved to {final_output}")

main()
