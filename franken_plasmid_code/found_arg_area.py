#this will take in the score of the filter file
# and the csv of args
# and tell you the area of the args that were found

import pandas as pd

def load_args_of_interest(file_with_args, skip_lines=8):
    with open(file_with_args, 'r') as f:
        lines = f.readlines()
    args = [line.strip().split(",") for line in lines[skip_lines:]]  # Get the ARG names from line 8 onwards
    #print(args)
    return args  # Using a set for fast lookup

def load_arg_locations(csv_file):
    df = pd.read_csv(csv_file)
    # Make sure positions are numeric
    df.iloc[:, 1] = pd.to_numeric(df.iloc[:, 1], errors='coerce')
    df.iloc[:, 2] = pd.to_numeric(df.iloc[:, 2], errors='coerce')

    # Build a dictionary: {arg_name: arg_length}
    arg_lengths = {}
    for idx, row in df.iterrows():
        arg_name = row.iloc[0]
        start = row.iloc[1]
        end = row.iloc[2]
        arg_lengths[arg_name] = end-start
    return arg_lengths

def calculate_found_bases(arg_lengths, args_of_interest):
    num_over_90 = 0
    output = []
    for arg, length in args_of_interest:
        arg_real_len = arg_lengths[arg]
        percent = float(length) / float(arg_real_len)
        if percent >=0.9:
            num_over_90 += 1
        output.append((arg,percent))


    return num_over_90, output

if __name__ == "__main__":
    args_file = "scored_franken_arg_blatdiver_overlap_div_114.tsv"
    locations_csv = "/usr1/gouallin/blat/franken_plasmid_tests/franken_arg/franken_arg.csv"

    args_of_interest = load_args_of_interest(args_file)
    arg_lengths = load_arg_locations(locations_csv)
    num_over_90, output = calculate_found_bases(arg_lengths, args_of_interest)

    print(f"Number of ARGs with ≥90% found: {num_over_90}")
    #print("\nARG recovery percentages:")
    #for arg, pct in output:
    #    print(f"{arg}: {pct:.2f}%")
