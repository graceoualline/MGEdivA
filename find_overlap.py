#this will only return regions that have been repeated, 
# and are different species

# will be O(N^2)
import csv
import sys
from filter_blat import *

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

    os.remove(q_seq)
    os.remove(ref_fasta_name)
    return distance

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None

#if it finds regions that overlap and map to the same species, it will combine them
def compress(input_file):
    #"Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7
    rows = set()

    with open(input_file, 'r') as f:
        #print("File", input_file, "opened")
        tsv = csv.reader(f, delimiter="\t")
        for line in tsv:
            if line[1].isnumeric():
                care = line[9:17] + line[21:] #q name to ref Name
                rows.add(tuple(care))

    #sort them by Qstart and Qend
    rows = sorted(list(rows), key=lambda x: (int(x[qs]), int(x[qe])))
    compressed = []
    #if there is no one, then return nothing
    if not rows:
        return []

    #Initialize first merged interval
    i = 0
    
    s_cur, e_cur, species_cur, row_cur = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], list(rows[i])
    while species_cur == "unclassified" and i < len(rows):
        # add those whose species is unknown
        compressed.append(tuple(row_cur))
        i += 1
        s_cur, e_cur, species_cur, row_cur = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], list(rows[i])
        
    
    for row in rows[i:]:
        s2, e2, species2, row2 = int(row[qs]), int(row[qe]), row[rsp], list(row)
        # if this new species is unclassified, just add it and move on
        if species2 == "unclassified":
            compressed.append(tuple(row2))
            continue
        
        if species2 == species_cur and overlaps(s_cur, e_cur, s2, e2):
            #then update the end
            e_cur = max(e_cur, e2)
        #if it does not overlap or does not share the same species, start the next merge
        else:
            #update this row to be the new values of the new end
            row_cur[qe] = str(e_cur)
            compressed.append(tuple(row_cur))
            #start a new current fam
            s_cur, e_cur, species_cur, row_cur = int(row[qs]), int(row[qe]), row[rsp], list(row)
        
    # add the final interval
    row_cur[qe] = str(e_cur)
    compressed.append(tuple(row_cur))
    #print(compressed)
    return compressed

def build_overlap_row(row1, row2, new_start, new_end, ani_value):
    new_row = []
    for h in range(len(row1)):
        if h in [0, 1, 9]:  # qname, rname, q_species
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
    qs = 2 #q start is 2
    qe = 3 # qend is 3
    rsp = 10 #ref species is 7
    new_rows = set()
    rows = sorted(rows, key=lambda x: (int(x[qs]), int(x[qe])))
    used = set()
    for i in range(len(rows)):
        if rows[i] in used:
            continue
        s1, e1, species1, row1 = int(rows[i][qs]), int(rows[i][qe]), rows[i][rsp], rows[i]
        #find the row with the most overlap that does not share a species
        for j in range(i+1, len(rows)):
            #print("j", j)
            if rows[j] in used:
                continue

            s2, e2, species2, row2 = int(rows[j][qs]), int(rows[j][qe]), rows[j][rsp], rows[j]
            new_start, new_end = max(s1, s2), min(e1, e2)

            if overlaps(s1, e1, s2, e2):
                if "unclassified" in [species1, species2]:
                    id1 = row1[4]
                    id2 = row2[4]
                    ani = find_ani_overlap(id1, id2, blat_db, gtdb_index)
                    if type(ani) != str and ani < 95:
                        new_row = build_overlap_row(row1, row2, new_start, new_end, ani)
                        new_rows.add(new_row)
                        used.add(rows[i])
                        used.add(rows[j])
                        break
                #if it overlaps and the species are not equal
                elif species1 != species2 and "unclassified" not in [species1, species2]:
                    #There was no ani between species considered since they had dif names
                    new_row = build_overlap_row(row1, row2, new_start, new_end, "NA")
                    new_rows.add(new_row)
                    used.add(rows[i])
                    used.add(rows[j])
                    break

    with open(output_file, "w") as out:
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tT size\tT start\tT end\tPercent Identity\tQuery Species\tReference Species\tDivergence Time\tANI bt seqs(if div=unk)\tANI bt ref seqs(if species unk)\n")
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

    rows = compress(input_file)
    find_overlap(rows, output, blat_db, gtdb_index)
    print("Filtered file written to", output)