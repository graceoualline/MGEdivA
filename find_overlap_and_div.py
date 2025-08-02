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
from find_overlap import *

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None

#if it finds regions that overlap and map to the same species, it will combine them

def get_div_alt(path1, path2, tree):
    #print(sp1, sp2)
    #path1 = get_path(sp1, tree)
    #print(sp1, sp2)
    if path1 == "NA" or path2 == "NA":
        return "unk:unable_to_find_ref_species_in_tree"

    distance = path1.get_distance(path2)
    
    return distance


def find_overlap_and_div(rows, output_file, tree, blat_db, kraken):
    #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7
    tname = 4
    new_rows = set()
    rows = sorted(list(rows), key=lambda x:int(x[qs]))
    used = set()
    species_path_cache = {}
    species_path_cache["unclassified"] = "NA"

    i = 0
    while i < len(rows):
        if rows[i] in used:
            i += 1
            continue

        s1, e1, species1, row1 = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], rows[i]
        # Cache species path
        if species1 not in species_path_cache:
            #if species1 == "unclassified_Arthrobacter": print("overlap_div", output_file)
            species_path_cache[species1] = get_path(species1, tree)
        
        j = i+1
        #find the row with the most overlap that does not share a species
        merged = False

        while j < len(rows):
            s2, e2, species2, row2 = int(rows[j][qs]), int(rows[j][qe]), rows[j][rsp], rows[j]
            if row2 in used:
                j += 1
                continue
            # if we completely surpassed where we can overlap, skip
            if s2 > e1:
                break
            #print(i, j)
            # if there is overlap, and the species or different or we dont know
            if overlaps(s1, e1, s2, e2) != None and (species1 != species2 or ("unclassified" in [species1, species2])):
                #check if species are different enough first
                if species2 not in species_path_cache:
                    species_path_cache[species2] = get_path(species1, tree)

                
                div = get_div_alt(species_path_cache[species1], species_path_cache[species2], tree)
                ani = "NA"
                if type(div) == str:
                    #get the ani instead
                    id1 = row1[tname]
                    id2 = row2[tname]
                    ani = find_ani_overlap(id1, id2, blat_db, kraken) #make find_ani later
                    #print("ani calculated:", ani)
                   
                if (type(div) != str and div >= 1) or (type(ani) == int and ani < 95):
                    #print("merging rows")
                    new_row = []
                    new_start = max(s1, s2)
                    new_end = min(e1, e2)
                    for h in range(len(row1)):
                        if h in [0, 1, rsp]: #0 and 1 are q name and length, and q species
                            new_row.append(row1[h])
                        elif h == 2: new_row.append(str(new_start)) #put new start or end
                        elif h==3: new_row.append(str(new_end)) # put in new end
                        else:
                            new_row.append(f"{row1[h]},{row2[h]}")
                    new_row.append(f"{div}")
                    new_row.append(f"{ani}")
                    #print(i, j, "merged")
                    used.add(rows[i])
                    used.add(rows[j])
                    new_rows.add("\t".join(new_row))
                    merged = True
                    break
            j+=1
        #this means row[i] had no takers, just keep looking
        if not merged:
            i += 1

    with open(output_file, "w") as out:
                 #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tDiv bt ref species\tAni bt ref species\n")
        #if we do have new rows, then write them
        if len(new_rows) > 0:
            for l in new_rows:
                out.write(l + "\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    #
    if len(sys.argv) != 6:
        print("Usage: python3 find_overlap_and_div.py <file of blatdiver output> <output file name> <tree> <blat_db> <index of species>")
        sys.exit(1)
    input_file = sys.argv[1]
    output = sys.argv[2]
    tree = Tree(sys.argv[3])
    blat_db = sys.argv[4]
    kraken = load_hash_table(sys.argv[5])

    print("Compressing overlap")
    rows = compress(input_file)
    print("Rows before:", len(rows))
    find_overlap_and_div(rows, output, tree, blat_db, kraken)
    print("Filtered file written to", output)