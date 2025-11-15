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

#if it finds regions that overlap and map to the same species, it will combine them

def find_overlap_and_div(rows, output_file, tree, blat_db, index):
    div_cache = dict()
    path_cache = dict()
    ani_cache = dict() # (id1, id2) : ani
    #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7
    tname = 4
    new_rows = set()
    rows = sorted(list(rows), key=lambda x:int(x[qs]))
    used = set()

    i = 0

    start_positions = [int(row[1][qs]) for row in rows]

    for i, row1 in enumerate(rows):
        if i in used:
            continue

        s1, e1, species1 = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp]
        
        if i + 1 < len(rows) and  e1 < int(rows[i+1][qs]):
            continue
        if species1 in path_cache: path1 = path_cache[species1]
        else: path1 = get_path_from_name(species1, tree)
        
        for j in range(i+1, len(rows)):
            #print(i, j)
            if j in used:
                continue
            s2, e2, species2, row2 = int(rows[j][qs]), int(rows[j][qe]), rows[j][rsp], rows[j]
            # if we completely surpassed where we can overlap, skip
            if s2 > e1:
                break

            if overlaps(s1, e1, s2, e2) == None:
                continue

            # if they are both equal and known
            if species1 == species2 and "unclassified" not in [species1, species2]:
                continue
            div = check_div_cache(species1, species2, div_cache)
            if div == None:
                if species2 in path_cache: path2 = path_cache[species2]
                else: path2 = get_path_from_name(species2, tree)
                # Get the divergence time between query species and reference species
                div = get_div(path1, path2, tree)
                div_cache[(species1, species2)] = div
            
            ani = "NA"

            # if div is not a number, then it is unknown
            if isinstance(div, str):
                #get the ani instead
                id1 = row1[tname]
                id2 = row2[tname]
                if (id1, id2) in ani_cache: ani = ani_cache[(id1, id2)]
                elif (id2, id1) in ani_cache: ani = ani_cache[(id2, id1)]
                else: 
                    ani = find_ani_overlap(id1, id2, blat_db, index)
                    ani_cache[(id1, id2)] = ani

                
            if (not isinstance(div, str) and div >= 1) or (isinstance(ani, int) and ani < 95):
                #print("merging rows")
                new_row = []
                new_start = max(s1, s2)
                new_end = min(e1, e2)
                for h in range(len(row1)):
                    if h in [0, 1, 9]: #0 and 1 are q name and length, and q species
                        new_row.append(row1[h])
                    elif h == 2: new_row.append(str(new_start)) #put new start or end
                    elif h==3: new_row.append(str(new_end)) # put in new end
                    else:
                        new_row.append(f"{row1[h]},{row2[h]}")
                new_row.extend([str(div), str(ani)])
                #print(i, j, "merged")
                used.add(rows[i])
                used.add(rows[j])
                new_rows.add("\t".join(new_row))
                break
            
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
    index = load_hash_table(sys.argv[5])

    #print("Compressing")
    rows = compress(input_file)
    #print("Rows before:", len(rows))
    find_overlap_and_div(rows, output, tree, blat_db, index)
    print("Filtered file written to", output)