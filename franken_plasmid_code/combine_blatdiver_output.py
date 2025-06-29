# this will combine all of the blatdiver output files into one big one so that we can do analysis on it
# should only be used at the end of everything so that there is less computation since it is N^2 to filter
# is only for franken arg analysis

import pandas as pd
import glob
import sys
import os

def adjust_and_merge_tsvs(chunk_dir, chunk_size, output_file, to_combine):
    #ex to_combine = *blatdiver_output.tsv
    all_files = sorted(glob.glob(os.path.join(chunk_dir, "*" + to_combine)))
    print(chunk_dir)
    print(to_combine)
    dfs = []
    total_len = 0
    for filename in all_files:
        print(filename)
        file = filename.split("/")[-1]
        idx = int(file.split("_")[1])
        # Calculate chunk offset
        offset = int(idx * chunk_size)

        with open(filename, 'r') as f:
            first_line = f.readline()
            has_header = not first_line.strip().split('\t')[0].isdigit()

        if has_header:
            df = pd.read_csv(filename, sep='\t')
        else:
            df = pd.read_csv(filename, sep='\t', header=None)
            
        total_len += len(df)
        
        if df.empty:
            continue  # Skip empty chunks

        # Assuming start is in column 2 and end is in column 3 (0-based)
        #change these depending on what file it is
        
        if df.shape[1] > 13:
            start_index = 11
            end_index = 12
        else:
            start_index = 2# 11 for full size
            end_index = 3
        #skip if its the first one, dont bother
        if offset > 0:
            df.iloc[:, start_index] = pd.to_numeric(df.iloc[:, start_index], errors = 'coerce') + offset  # Adjust start
            df.iloc[:, end_index] = pd.to_numeric(df.iloc[:, end_index], errors = 'coerce') + offset  # Adjust end

        dfs.append(df)
        
        

    # Merge all
    final_df = pd.concat(dfs, ignore_index=True)

    def is_misaligned(row):
            # First two columns are blank but middle columns are not
            return pd.isna(row[0]) and pd.isna(row[1]) and not pd.isna(row[10])

    bad_rows = final_df[final_df.apply(is_misaligned, axis=1)]
    print("Found misaligned rows:")
    print(bad_rows.to_string(index=False))

    # Save
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f"Combined file written to {output_file}")
    print("Total length is", total_len)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 combine_blatdiver_output.py <blatdiver_output_directory> <chunk_size> <which files to combine ex. *blatdiver_output.tsv> <output_file>")
        sys.exit(1)

    chunk_dir = sys.argv[1]
    chunk_size = int(sys.argv[2])
    to_combine = sys.argv[3]
    output_file = sys.argv[4]

    adjust_and_merge_tsvs(chunk_dir, chunk_size, output_file, to_combine)
