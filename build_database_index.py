# this should intake the folder that contains a bunch of split of 2bit files
# and it will return a .pkl file which is a hash table that essential works like
# index[sequence_id] = (sequence_species, sequence_location_in_the_database)
import os
import subprocess
import tempfile
import pickle
from extract_species_from_kraken import *

#first step, get the locations of all sequence ids
def get_seq_id_loc(two_bit_dir, input_files):
    seq_id_loc = set()
    
    print(f"Processing {len(input_files)} 2bit files...")
    for i, inputf in enumerate(input_files, 1):
        two_bit_path = os.path.join(two_bit_dir, inputf)

        print(f"Progress: {i}/{len(input_files)} - Processing {inputf}")

        #this will make a temp file
        with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
            temp_path = temp_file.name

        # Run twoBitInfo to get sequence info
        try:
            result = subprocess.run(["twoBitInfo", two_bit_path, temp_path],
                                    capture_output = True,
                                    text = True,
                                    check = True)
        
            # Read the sequence IDs and write with file name
            with open(temp_path, "r") as f:
                for line in f:
                    seq_id = line.split()[0]
                    seq_id_loc.add((seq_id, inputf))
        except subprocess.CalledProcessError as e:
            print(f"ERROR processing {inputf}: {e}")
        finally:
            os.remove(temp_path)
    print(f"Found {len(seq_id_loc)} total sequences")
    return seq_id_loc

def combine_seq_spec_loc(species_file, seq_id_loc, two_bit_dir):
    # Step 1: Load species data into a dictionary
    seq_to_species = {}
    if species_file and os.path.exists(species_file):
        print(f"Loading existing species data from {species_file}...")
        with open(species_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    seq_id, species = parts[0], parts[1]
                    seq_to_species[seq_id] = species
        print(f"Loaded {len(seq_to_species)} existing species annotations")
    else:
        raise Exception("No species file provided or file doesn't exist")

    #Step 2: Process each sequence
    result_dict = {}
    
    for i, seq_loc in enumerate(seq_id_loc, 1):
        #make sure there are two entries
        if len(seq_loc) >= 2:
            seq_id, file_location = seq_loc[0], seq_loc[1]

            # Progress indicator
            if i % 10000 == 0 or i == len(seq_id_loc):
                print(f"Progress: {i}/{len(seq_id_loc)} sequences processed")

            # if the sequence does not have a species, go find it with kraken
            species = seq_to_species.get(seq_id, "unclassified")
            result_dict[seq_id] = (species, file_location)

    return result_dict

def load_hash_table(index_file):
    """
    Loads the precomputed hash table from a binary file.
    :param index_file: Path to the stored hash table
    :return: Dictionary {seq_id: species}
    """
    with open(index_file, "rb") as f:
        return pickle.load(f)

def lookup_species(hash_table, seq_id):
    """
    Looks up the species name for a given seq_id using the hash table.
    Returns "unclassified" if not found.
    """
    tuple = hash_table.get(seq_id, "unclassified")
    return tuple[0]

def lookup_location(hash_table, seq_id):
    """
    Looks up the species name for a given seq_id using the hash table.
    Returns "unclassified" if not found.
    """
    tuple = hash_table.get(seq_id, "unclassified")
    return tuple[1]

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 build_database_index.py <2bit_blat_db_dir> <output_file_name> <species_file>")
        print("  2bit_blat_db_dir: Directory containing .2bit files")
        print("  output_file_name: Output .pkl file name")
        print("  species_file: Existing species annotation file")
        sys.exit(1)
    two_bit_dir = sys.argv[1]
    final_output = sys.argv[2]
    species_file = sys.argv[3]

    if not os.path.exists(two_bit_dir):
        print(f"ERROR: 2bit directory {two_bit_dir} does not exist")
        sys.exit(1)

    input_files = [f for f in os.listdir(two_bit_dir) 
                 if f.endswith('.2bit') and os.path.isfile(os.path.join(two_bit_dir, f))]

    if not input_files:
        print(f"ERROR: No .2bit files found in {two_bit_dir}")
        sys.exit(1)

    #first step, get the locations of all sequence ids
    #used to be get_seq_ids.py
    print("\n" + "="*60)
    print("STEP 1: Extracting sequence IDs from 2bit files")
    print("="*60)
    seq_id_loc = get_seq_id_loc(two_bit_dir, input_files)
    # next, we use a file that connects seq_id and species determined by kraken
    # and add that to the index # use to be combine_species_location.py
    print("\n" + "="*60) 
    print("STEP 2: Combining with species data and Kraken classification")
    print("="*60)
    seq_species_loc = combine_seq_spec_loc(species_file, seq_id_loc, two_bit_dir)

    # # now make a loadable index from this new table of id, species, location
    print("\n" + "="*60)
    print("STEP 3: Saving database index")
    print("="*60)
    with open(final_output, "wb") as f:
        pickle.dump(seq_species_loc, f)

    print(f"Hash index saved to {final_output}")
    print("="*60)
    print("DATABASE INDEX CREATION COMPLETED!")
    print("="*60)
