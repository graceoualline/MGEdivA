#this will only return regions that have been repeated, 
# and are different species

# will be O(N^2)
import csv
import sys
from find_overlap_and_div import *

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None

#if it finds regions that overlap and map to the same species, it will combine them
def compress(input_file, output_file):
    rows = set()

    with open(input_file, 'r') as f:
        print("File", input_file, "opened")
        tsv = csv.reader(f, delimiter="\t")
        for line in tsv:
            if line[1].isnumeric():
                care = line[9:15] + line[21:] #q name to ref Name
                rows.add(tuple(care))

    rows = sorted(list(rows), key=lambda x: (int(x[2]), int(x[3])))
    compressed = []
    #if there is no one, then return nothing
    if not rows:
        return []

    #Initialize first merged interval
    i = 0
    s_cur, e_cur, species_cur, row_cur = int(rows[i][2]), int(rows[i][3]), rows[i][7], list(rows[i])
    while species_cur == "unclassified" and i < len(rows):
        # add those whose species is unknown
        compressed.append(tuple(row_cur))
        i += 1
        s_cur, e_cur, species_cur, row_cur = int(rows[i][2]), int(rows[i][3]), rows[i][7], list(rows[i])
        
    
    for row in rows[i:]:
        s2, e2, species2, row2 = int(row[2]), int(row[3]), row[7], list(row)
        # if this new species is unclassified, just add it and move on
        if species2 == "unclassified":
            compressed.append(tuple(row2))
            continue
        
        #if it does not overlap or does not share the same species, start the next merge
        if species2 != species_cur or not overlaps(s_cur, e_cur, s2, e2):
            #update this row to be the new values of the new end
            row_cur[3] = str(e_cur)
            compressed.append(tuple(row_cur))
            #start a new current fam
            s_cur, e_cur, species_cur, row_cur = int(row[2]), int(row[3]), row[7], list(row)
        elif species2 == species_cur and overlaps(s_cur, e_cur, s2, e2):
            #then update the end
            e_cur = max(e_cur, e2)
        else:
            print("Error, weird third case")
            assert(False)

    # add the final interval
    row_cur[3] = str(e_cur)
    compressed.append(tuple(row_cur))
    #print(compressed)
    compressed_output = f"compressed_{output_file}"
    with open(compressed_output, "w") as out:
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tTsize\tQuery_Species\tReference_Species\tDivergence_Time\tANI_of_whole_seqs(only_if_div=unk)\n")
        for l in compressed:
            line = "\t".join(l)
            out.write(line + "\n")
    return compressed

def build_overlap_row(row1, row2, new_start, new_end, ani_value):
    new_row = []
    for h in range(len(row1)):
        if h in [0, 1, 6]:  # qname, rname, q_species
            new_row.append(row1[h])
        elif h == 2:
            new_row.append(str(new_start))
        elif h == 3:
            new_row.append(str(new_end))
        else:
            new_row.append(f"{row1[h]},{row2[h]}")
    new_row.append(str(ani_value))
    return "\t".join(new_row)

def find_overlap(rows, output_file, blat_db, gtdb_index):
    new_rows = set()
    rows = sorted(rows, key=lambda x: (int(x[2]), int(x[3])))
    used = set()
    for i in range(len(rows)):
        if rows[i] in used:
            continue
        s1, e1, species1, row1 = int(rows[i][2]), int(rows[i][3]), rows[i][7], rows[i]
        #find the row with the most overlap that does not share a species
        for j in range(i+1, len(rows)):
            #print("j", j)
            if rows[j] in used:
                continue

            s2, e2, species2, row2 = int(rows[j][2]), int(rows[j][3]), rows[j][7], rows[j]
            new_start, new_end = max(s1, s2), min(e1, e2)

            if overlaps(s1, e1, s2, e2):
                if "unclassified" in [species1, species2]:
                    id1 = row1[4]
                    id2 = row2[4]
                    # goes to find_overlap_and_div.py
                    ani = find_ani_overlap(id1, id2, blat_db, gtdb_index)
                    if ani != None and ani < 95:
                        new_row = build_overlap_row(row1, row2, new_start, new_end, ani)
                        new_rows.add(new_row)
                        used.add(rows[i])
                        used.add(rows[j])
                        break
                #if it overlaps and the species are not equal
                elif species1 != species2 and "unclassified" not in [species1, species2]:
                    #There was no ani between species considered since they had dif names
                    new_row = new_row = build_overlap_row(row1, row2, new_start, new_end, "NA")
                    new_rows.add(new_row)
                    used.add(rows[i])
                    used.add(rows[j])
                    break

    with open(output_file, "w") as out:
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tTsize\tQuery_Species\tReference_Species\tDivergence_Time\tANI_of_whole_seqs(only_if_div=unk)\tANI_bt_r_seqs(if_species_unk)\n")
        for l in new_rows:
            out.write(l + "\n")

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    #
    if len(sys.argv) != 5:
        print("Usage: python3 find_overlap.py <file of blatdiver output> <output file name> <blat_db> <gtdb_species_index>")
        sys.exit(1)
    input_file = sys.argv[1]
    output = sys.argv[2]
    blat_db = sys.argv[3]
    gtdb_index = load_hash_table(sys.argv[4])

    rows = compress(input_file, output)
    find_overlap(rows, output, blat_db, gtdb_index)
    print("Filtered file written to", output)