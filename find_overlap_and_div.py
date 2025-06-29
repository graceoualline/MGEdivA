#same as find overlap, but we also filter for if the two lines we are compressing 
# map to two different speceis that are divergently different

#this will only return regions that have been repeated, 
# and are different species

# will be O(N^2)
import csv
import sys
from filter_blat import *
from ete3 import Tree
import os

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None


def find_ani_overlap(q_id, ref_id, blat_db, kraken):
    """
    Finds the ANI between a query sequence and a reference sequence stored in a .2bit file.

    Parameters:
    - q_seq (str): Path to the query sequence (FASTA format).
    - ref_id (str): The reference sequence ID to extract.
    - blat_db (str): Path to the folder containing .2bit files.

    Returns:
    - float: ANI value (or None if reference is not found).
    """
    q_seq = extract_2bit_fasta(q_id, kraken, blat_db)
    ref_fasta_name = extract_2bit_fasta(ref_id, kraken, blat_db)
        
    distance = calculate_distance(q_seq, ref_fasta_name)
    
    os.remove(ref_fasta_name)
    return distance

#if it finds regions that overlap and map to the same species, it will combine them
def compress(input_file):
    rows = set()

    with open(input_file, 'r') as f:
        print("File", input_file, "opened")
        tsv = csv.reader(f, delimiter="\t")
        for line in tsv:
            #print(line)
            if line[1].isnumeric():
                care = line[9:15] + line[21:] #q name to ref Name
                rows.add(tuple(care))
            #print("care", care)
    print("Number of lines before:", len(rows))

    counter = 0
    rows = list(rows)
    rows.sort()
    new_rows = set()
    used = set()
    i = 0
    
    while i < len(rows):
        #print("i", i)
        if rows[i] in used:
            i += 1
            continue
        #print(rows[i])
        s1, e1, species1, row1 = int(rows[i][2]), int(rows[i][3]), rows[i][7], rows[i]
        #print(row1)
        if species1 == "unclassified":
            i += 1
            continue
        j = i + 1
        merged = False
        while j < len(rows):
            #print("j", j)
            if rows[j] in used:
                j+=1
                continue

            s2, e2, species2, row2 = int(rows[j][2]), int(rows[j][3]), rows[j][7], rows[j]
            if species2 == "unclassified":
                j += 1
                continue
            if overlaps(s1, e1, s2, e2) != None and species1 == species2:
                new_row = ""
                new_start = min(s1, s2)
                new_end = max(e1, e2)
                for h in range(len(row1)):
                    if h in [0, 1, 4, 5, 6, 7, 8, 9]: #these should not change as we compress
                        new_row += f"{row1[h]}\t"
                    elif h == 2: #put new start or end
                        new_row += f"{new_start}\t"
                    elif h==3: # put in new end
                        new_row += f"{new_end}\t"
                new_row = new_row[:-1] #remove the last \t
                new_row = new_row.split("\t")
                counter += 1
                used.add(rows[i]) #remove this line since we already merged it
                used.add(rows[j]) #remove this line since we already merged it
                new_rows.add(tuple(new_row))
                merged = True
                break
            j += 1
        #if there was no one to merge with, skip it
        if merged == False:
            i += 1
    #print("number of lines after compression:", len(new_rows))
    #with open("rand97_hidden_plasmid_compressed.tsv", "w") as out:
    #    out.write("Q name\tQ size\tQ start\tQ end\tT name\tTsize\tQuery_Species\tReference_Species\tDivergence_Time\tANI_of_whole_seqs(only_if_div=unk)\n")
    #    for l in new_rows:
    #        row = "\t".join(l)
    #        out.write(row + "\n")

    return list(new_rows)

def get_div_alt(path1, path2, tree):
    #print(sp1, sp2)
    #path1 = get_path(sp1, tree)
    #print(sp1, sp2)
    if path1 == None or path2 == None:
        return None

    distance = path1.get_distance(path2)
    
    return distance


def find_overlap_and_div(rows, output_file, tree, blat_db, kraken):
    counter = 0
    i = 0
    new_rows = set()
    rows = sorted(list(rows), key=lambda x:int(x[2]))
    rows.sort()
    used = set()
    #print("ROW 399", rows[399])
    #assert(False)
    species_path_cache = {}

    row_index = 0
    while i < len(rows):
        if rows[i] in used:
            i += 1
            continue
        #print("i", i)
        s1, e1, species1, row1 = int(rows[i][2]), int(rows[i][3]), rows[i][7], rows[i]
        # Cache species path
        if species1 not in species_path_cache:
            species_path_cache[species1] = get_path(species1, tree)
        
        #while we the end doesnt surpass the start
        while row_index < len(rows) and int(rows[i][3]) < s1:
            row_index += 1
        j = row_index

        #find the row with the most overlap that does not share a species
        merged = False
        # as long as the start of an arg is less than the end:
        while j < len(rows) and int(rows[j][2]) <= e1:
            #print("j", j)
            if rows[j] in used:
                j += 1 
                continue
            s2, e2, species2, row2 = int(rows[j][2]), int(rows[j][3]), rows[j][7], rows[j]

            if overlaps(s1, e1, s2, e2) != None and (species1 != species2 or (species1 == "unclassified" or species2 == "unclassified")):
                #check if species are different enough first
                div = None
                ani = None
                if species2 not in species_path_cache:
                    species_path_cache[species2] = get_path(species1, tree)
                div = get_div_alt(species_path_cache[species1], species_path_cache[species2], tree)
                #print(f"Divergence {species1} {species2} {div}")
                if div == None:
                    #get the ani instead
                    id1 = row1[4]
                    id2 = row2[4]
                    ani = find_ani_overlap(id1, id2, blat_db, kraken) #make find_ani later
                    #print("ani calculated:", ani)
                   
                if (div != None and div >= 1) or (ani != None and ani < 95):
                    #print("merging rows")
                    new_row = ""
                    new_start = max(s1, s2)
                    new_end = min(e1, e2)
                    for h in range(len(row1)):
                        if h in [0, 1, 6]: #0 and 1 are q name and length, and 5 is q species
                            new_row += f"{row1[h]}\t"
                        elif h == 2: #put new start or end
                            new_row += f"{new_start}\t"
                        elif h==3: # put in new end
                            new_row += f"{new_end}\t"
                        else:
                            new_row = new_row + f"{row1[h]},{row2[h]}\t"
                    new_row = new_row + f"{div}\t{ani}"
                    counter += 1
                    used.add(rows[i])
                    used.add(rows[j])
                    rows.remove(rows[j]) #remove this line since we already merged it
                    new_rows.add(new_row)
                    merged = True
                    break
            j+=1
        if merged == False:
            i += 1

    with open(output_file, "w") as out:
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tTsize\tQuery_Species\tReference_Species\tDivergence_Time_to_Q\tANI_of_whole_seqs(only_if_div=unk)\tDiv_bt_species\tAni_bt_species\n")
        #if we do have new rows, then write them
        if len(new_rows) > 0:
            for l in new_rows:
                out.write(l + "\n")
    print("number of lines before filter:", len(rows))
    print("number of lines after overlap:", len(new_rows))

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    #
    if len(sys.argv) != 6:
        print("Usage: python3 find_overlap_and_div.py <file of blatdiver output> <output file name> <tree> <blat_db> <kraken>")
        sys.exit(1)
    input = sys.argv[1]
    output = sys.argv[2]
    tree = Tree(sys.argv[3])
    blat_db = sys.argv[4]
    kraken = load_hash_table(sys.argv[5])

    rows = compress(input)
    find_overlap_and_div(rows, output, tree, blat_db, kraken)
    print("Filtered file written to", output)