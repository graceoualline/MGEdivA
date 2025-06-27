# will just do the following:
#If Q_gap_bases > matches: dont include it
#if this works well I will just fold this in the filter_blat code

import sys
import csv

def remove_large_gaps(input_file, output_file):
    with open(input_file, 'r') as f, open(output_file, 'w', newline='') as out:
        print("File", input_file, "opened")
        tsv_reader = csv.reader(f, delimiter="\t")
        tsv_writer = csv.writer(out, delimiter="\t")
        count_old = 0
        count_new = 0

        # Read and write the header
        header = next(tsv_reader)
        tsv_writer.writerow(header)
        
        for line in tsv_reader:
            count_old += 1
            #make sure is not a header
            if line[1].isnumeric():
                # if number of matches >= number of gaps in the Query: keep the line
                if int(line[0]) >= int(line[5]):
                    tsv_writer.writerow(line)
                    count_new += 1
        print("Number of lines before:", count_old)
    print("Number of lines after:", count_new)

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    # if q_species = unk, then run the get query species
    if len(sys.argv) != 3:
        print("Usage: python3 remove_large_gaps.py <blatdiver_output.tsv> <output file name>")
        sys.exit(1)

    # Usage example:
    input_file = sys.argv[1]
    output_name = sys.argv[2]
    
    remove_large_gaps(input_file, output_name)
 