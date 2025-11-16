# this should intake the folder that contains a bunch of split of 2bit files
# and it will return a .pkl file which is a hash table that essential works like
# index[sequence_id] = (sequence_species, sequence_location_in_the_database, name_in_tree)
from ete3 import Tree
import os
import subprocess
import tempfile
import pickle
from extract_species_from_kraken import *
import pandas as pd
import numpy as np

def get_one(path_iter):
    paths = list(path_iter)
    if len(paths) != 1:
        raise ValueError(f"Expected 1 file, found {len(paths)}: {paths}. You do not have the proper files needed for the divergence tree.")
    return str(paths[0])

def load_H5(file):
    store = pd.HDFStore(file)
    df = store["table"]
    return df

class Divergence_Tree_Preprocessed():
    """
    Class for finding divergence along a tree that has been preprocessed
    """
    def __init__(self, tree_dir):
        tree_path = Path(tree_dir)
        nwk_file   = get_one(tree_path.glob("*.nwk"))
        plk_file   = get_one(tree_path.glob("*.pkl"))
        index_file = get_one(tree_path.glob("*.index.h5"))
        mins_file  = get_one(tree_path.glob("*.mins.h5"))
        tour_file  = get_one(tree_path.glob("*.tour.h5"))

        print("Loading Divergence Tree")
        self.tree = Tree(nwk_file)
        self.species_dist = load_hash_table(plk_file)
        self.index = load_H5(index_file)
        self.mins = load_H5(mins_file)
        self.tour = load_H5(tour_file)[0].tolist()
        print("Divergence Tree Loaded!")
        
    def get_min_index(self, i,j):
        """Range min query on L[i:j]"""
        m = int(np.log2(j-i))
        minimum = min(self.mins.iat[m, i], self.mins.iat[m, j-2**m])
        min_index = self.index.iat[m, i] if self.mins.iat[m, i] == minimum else self.index.iat[m, j-2**m]
        return min_index

    def divergence(self, a, b):
        """Finds divergence between nodes a and b"""
        if a == "NA" or b == "NA":
            return "unk:unable_to_find_ref_species_in_tree"
        a_tour_step = self.species_dist[a][1]
        b_tour_step = self.species_dist[b][1]
        a_dist = self.species_dist[a][0]
        b_dist = self.species_dist[b][0]
        lca_step = self.get_min_index(min(a_tour_step, b_tour_step), max(a_tour_step, b_tour_step)+1)
        lca = self.tour[lca_step]
        lca_dist = self.species_dist[lca][0]
        return (a_dist + b_dist - 2*lca_dist)/2

    def get_path(self, sp):
        if sp == "unclassified":
            return "NA"

        sp = "_".join(sp.split(" "))
        sp = "'" + sp + "'"

        try:
            path = self.tree & sp
            return path
        except:
            first = sp.split("_")
            try:
                new = first[0] + "_" + first[1] + "'"
                path = self.tree & new
                return path
            except:
                try: 
                    new = first[0] + "'"
                    path = self.tree & new
                    return path
                except:
                    #if all else fails, find the first instance of the name within another name
                    for node in self.tree.traverse("preorder"):
                        if first[0] in node.name:
                            return node
                    return "NA"

    # return the actual node name that is within the tree
    def get_path_identifier(self, node):
        if node == "NA":
            return "NA"
        path_names = []
        current = node
        while current is not None:
            path_names.append(current.name if current.name else "unnamed")
            current = current.up
        return path_names[0]

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

def combine_seq_spec_loc(species_file, seq_id_loc, two_bit_dir, tree):
    # Step 1: Load species data into a dictionary
    seq_to_species = {}
    all_species = set()
    if species_file and os.path.exists(species_file):
        print(f"Loading existing species data from {species_file}...")
        with open(species_file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 2:
                    seq_id, species = parts[0], parts[1]
                    all_species.add(species)
                    seq_to_species[seq_id] = species
        print(f"Loaded {len(seq_to_species)} existing species annotations")
    else:
        raise Exception("No species file provided or file doesn't exist")
    ##Debug
    # seq_to_species_test = [('JAHAAN010000031.1', 'unclassified'), ('QHWJ01000176.1', 'Rhizobiaceae'), ('WQSI01000087.1', 'unclassified'), ('NZ_QJUG01000292.1', 'Nonomuraea gerenzanensis'), ('JAJXTT010000076.1', 'Paracoccus versutus'), ('JAADEE010000229.1', 'Acinetobacter baumannii')]
    # all_species = [i for _, i in seq_to_species_test]
    # print(seq_to_species_test)
    # print(all_species)
    #Step 2: Create a dictionary of all species, and where they are in the tree
    print("Obtaining species path in the tree")
    print(f"There are {len(all_species)} species")
    species_leaf_name = dict()
    for i, sp in enumerate(all_species):
        if i % 50 == 0:
            print(f"{i}/{len(all_species)} paths found.")
        path = tree.get_path(sp)
        leaf_name = tree.get_path_identifier(path)
        species_leaf_name[sp] = leaf_name
    # print("Species leaf name dictionary", species_leaf_name)
    print("Paths obtained. Creating hash table.")
    #Step 3: Process each sequence
    result_dict = {}
    # print("seq_id_loc:", seq_id_loc)
    for i, seq_loc in enumerate(seq_id_loc, 1):
        #make sure there are two entries
        if len(seq_loc) >= 2:
            seq_id, file_location = seq_loc[0], seq_loc[1]

            # Progress indicator
            if i % 100 == 0 or i == len(seq_id_loc):
                print(f"Progress: {i}/{len(seq_id_loc)} sequences processed")

            # if the sequence does not have a species, go find it with kraken
            species = seq_to_species.get(seq_id, "unclassified")
            #print("seq_id, seq_to_species[seq_id]", seq_id, species)
            result_dict[seq_id] = (species, file_location, species_leaf_name[species])
    # print("RESULT DIcT:", result_dict)
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

def lookup_tree_leaf_name(hash_table, seq_id):
    """
    Looks up the species name for a given seq_id using the hash table.
    Returns "unclassified" if not found.
    """
    tuple = hash_table.get(seq_id, "NA")
    return tuple[2]

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 build_database_index.py <2bit_blat_db_dir> <output_file_name> <species_file> <tree>")
        print("  2bit_blat_db_dir: Directory containing .2bit files")
        print("  output_file_name: Output .pkl file name")
        print("  species_file: Existing species annotation file")
        print("  tree: Time Tree of Life .nwk file")
        sys.exit(1)
    two_bit_dir = sys.argv[1]
    final_output = sys.argv[2]
    species_file = sys.argv[3]
    tree = Divergence_Tree_Preprocessed(sys.argv[4])

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
    #seq_id_loc = [('JAHAAN010000031.1', 'f1'), ('QHWJ01000176.1', 'f2'), ('WQSI01000087.1', 'f3'), ('NZ_QJUG01000292.1', 'f4'), ('JAJXTT010000076.1', 'f5'), ('JAADEE010000229.1', 'f6')]
    # next, we use a file that connects seq_id and species determined by kraken
    # and add that to the index # use to be combine_species_location.py
    print("\n" + "="*60) 
    print("STEP 2: Combining with species data and Kraken classification")
    print("="*60)
    seq_species_loc = combine_seq_spec_loc(species_file, seq_id_loc, two_bit_dir, tree)

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
