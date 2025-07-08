#code for making a hash table lookup for sequence ids so we can quickly find their species
# and for searching the hash table

import pickle

def build_and_store_hash_table(kraken_file, index_file):
    """
    Reads a file and stores the hash table (dictionary) to a binary file.
    :param kraken_file: Path to the Kraken file (tab-separated: seq_id -> species)
    :param index_file: Path to store the hash table
    """
    hash_table = {}

    with open(kraken_file, "r") as file:
        for line in file:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue  # Skip malformed lines
            
            seq_id, species, location = parts[0], parts[1], parts[2]
            hash_table[seq_id] = (species, location)  # Store in hash table

    # Save the hash table to a file using pickle
    with open(index_file, "wb") as f:
        pickle.dump(hash_table, f)

    print(f"Hash index saved to {index_file}")


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

# Example usage
# Example usage
#kraken_file = "/usr1/gouallin/blat/all_gtdb_seq_species_location_smart_split.txt"  # Input file
#index_file = "/usr1/shared/all_gtdb_seq_species_location_index_smart_split.pkl"  # Output file for the hash table
#build_and_store_hash_table(kraken_file, index_file)
'''#tests
hash_table = load_hash_table(index_file)  # Load the hash table once
print("loaded")
query_seq_id = "NZ_CP017181.1"
species = lookup_species(hash_table, query_seq_id)
print(f"Species for {query_seq_id}: {species}")
species = lookup_location(hash_table, query_seq_id)
print(f"location for {query_seq_id}: {species}")
'''
