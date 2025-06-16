#this will only return regions that have been repeated, 
# and are different species

# will be O(N^2)
import csv
import sys

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

import csv
import sys

#will keep Q name, Q start, Q end, T name(s), Q spec, T specie(s), divergence time(s), ani

def overlaps(start1, end1, start2, end2):
    if max(start1, start2) <= min(end1, end2):
        return min(end1, end2) - max(start1, start2)
    return None

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


def find_overlap(rows, output_file):
    counter = 0
    i = 0
    new_rows = set()
    rows = list(rows)
    rows.sort()
    used = set()
    while i < len(rows):
        if rows[i] in used:
            i += 1
            continue
        #print("i", i)
        s1, e1, species1, row1 = int(rows[i][2]), int(rows[i][3]), rows[i][7], rows[i]
        #FIX LATER I am ignoring unclassified species, should instead look at ANI
        if species1 == "unclassified":
            i += 1
            continue
        #find the row with the most overlap that does not share a species
        j = i + 1
        merged = False
        while j < len(rows):
            #print("j", j)
            if rows[j] in used:
                j += 1 
                continue

            s2, e2, species2, row2 = int(rows[j][2]), int(rows[j][3]), rows[j][7], rows[j]
            if species2 == "unclassified":
                j += 1
                continue
            if overlaps(s1, e1, s2, e2) != None and species1 == species2:
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
                new_row = new_row[:-1] #remove the last \t
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
        out.write("Q name\tQ size\tQ start\tQ end\tT name\tTsize\tQuery_Species\tReference_Species\tDivergence_Time\tANI_of_whole_seqs(only_if_div=unk)\n")
        for l in new_rows:
            out.write(l + "\n")
    print("number of lines after overlap:", len(new_rows))

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    #
    if len(sys.argv) != 3:
        print("Usage: python3 find_overlap.py <file of blatdiver output> <output file name>")
        sys.exit(1)
    input = sys.argv[1]
    output = sys.argv[2]
    rows = compress(input)
    find_overlap(rows, output)
    print("Filtered file written to", output)