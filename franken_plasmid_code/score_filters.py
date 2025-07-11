# this will intake the csv of arg locations, and .tsv of a filter, 
# and output the number of args found, area of args found, area of false pos
# and total area of args

import pandas as pd
import sys


#this should just return the Q start and Q end columns as its own df
def load_regions_from_tsv(tsv_file):
    #df = pd.read_csv(tsv_file, sep='\t', header = None)
    df = pd.read_csv(tsv_file, sep='\t')
    # Assuming the TSV has 'start' and 'end' in columns 2 and 3 (zero-indexed as 1 and 2)
    print(df.columns.tolist())
    #new_df = df[[11, 12]].copy()
    try:
        new_df = df[['Q start', 'Q end']].copy()
    except:
        new_df = df[['Q start', ' Q end']].copy()
    return new_df

# get rid of overlapping intervals so we aren't double counting the score
def merge_intervals(intervals_df):
    """Merge overlapping intervals."""
    
    if intervals_df.empty:
        return []
    
    # Convert DataFrame to list of tuples and sort
    intervals = list(zip(intervals_df.iloc[:, 0], intervals_df.iloc[:, 1]))
    sorted_intervals = sorted(intervals)
    
    merged = [sorted_intervals[0]]
    
    for current_start, current_end in sorted_intervals[1:]:
        last_start, last_end = merged[-1]
        
        if current_start <= last_end:  # Overlap
            merged[-1] = (last_start, max(last_end, current_end))
        else:
            merged.append((current_start, current_end))
    
    # Convert back to DataFrame
    return merged

def region_overlap(r1, r2):
    return max(0, min(float(r1[1]), float(r2[1])) - max(float(r1[0]), float(r2[0])))

def calculate_coverage(coor, arg_regions):
    found_args = {}
    total_overlap_area = 0
    total_found_area = 0
    
    #optimize this to only look at areas where it may overlap, do that later
    arg_index = 0
    # i is index, c is start and end
    for i, c in enumerate(coor): 
        #print(i)
        c_start, c_end = int(c[0]), int(c[1])
        total_found_area += c_end - c_start
        # Skip arg regions that are completely before this coordinate
        # as long as the end of the arg is less than the start of this region, skip
        while arg_index < len(arg_regions) and int(arg_regions.iloc[arg_index, 2]) < c_start:
            arg_index += 1
        j = arg_index
        # as long as the start of an arg is less than the end:
        while j < len(arg_regions) and int(arg_regions.iloc[j, 1]) <= c_end:
            start = arg_regions.iloc[j, 1]
            end = arg_regions.iloc[j, 2]
            overlap = region_overlap(c, (start, end))
            if overlap > 0:
                arg = arg_regions.iloc[j, 0]
                if arg in found_args: found_args[arg] += overlap
                else: found_args[arg] = overlap

                total_overlap_area += overlap
            j += 1
    return found_args, total_overlap_area, total_found_area - total_overlap_area

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python3 score_filters.py <arg csv> <output filtered tsv file> <output metrics file name>")
        sys.exit(1) 
    arg_csv = sys.argv[1]  # your CSV of inserted ARGs
    found_tsv = sys.argv[2]  # your software output
    output = sys.argv[3]

    df = pd.read_csv(arg_csv)
    #remove any args that exist after 11500000 bp because that is when it prematurely stopped
    #arg_regions = df[df.iloc[:, 1] <= 11500000].copy()  # REMOVE THIS IF WE END UP DOING MORE THAN 114 CHUNKS, filtered for args within the first 114 chunks
    #this is first column, the arg name, start, end
    arg_regions = df
    #print((f"Area of args that were inserted: {(arg_regions.iloc[:, 2] - arg_regions.iloc[:, 1]).sum()}"))


    # this is a df with the Q start and Q end sections of 
    found_regions = load_regions_from_tsv(found_tsv)
    print("loaded the tsv")

    #if there are any intervals found that overlap, this gets rid of them
    merged_found_regions = merge_intervals(found_regions)
    print("num regions found", len(merged_found_regions))

    # found args is a dictionary of which arg was found, and its amount of overlap
    # within_args is just the total area found within the args
    # outside_args is total area found outside args
    found_args, within_args, outside_args = calculate_coverage(merged_found_regions, arg_regions)


     #this is the total area of all args found
    area_of_all_args = (arg_regions.iloc[:, 2] - arg_regions.iloc[:, 1]).sum()

    #this is adding an area column for each arg
    arg_regions['length_bp'] = arg_regions.iloc[:, 2] - arg_regions.iloc[:, 1]

    #this is mapping the found arg dictionary to the names in the arg regions
    arg_regions['found_bp'] = arg_regions.iloc[:, 0].map(found_args)

    # get the percent found of that arg
    arg_regions['percent_area_found'] = arg_regions['found_bp'] / arg_regions['length_bp']
    #this is the total area we found
    found_area_of_args = arg_regions['found_bp'].sum()
    #this is the total area we should have found

    all_area_of_found_args = 0
    for arg in found_args:
        row = arg_regions[arg_regions.iloc[:, 0] == arg]
        if not row.empty:
            all_area_of_found_args += int((row.iloc[0, 2] - row.iloc[0, 1]))

    #get the number of args where >= 90% of their area was found
    num_over_90 = (arg_regions['percent_area_found'] >= 0.9).sum()

    with open(output, 'w') as output:
        output.write(f"Number of ARGs found: {len(found_args)}\n")
        output.write(f"Number of ARGs that were inserted: {len(arg_regions)}\n")
        #get the area of the args that were inserted
        output.write(f"percentage of args found: {len(found_args)/len(arg_regions)}\n")
        #REPLACE WITH TOTAL AREA OF HOST GENOME
        output.write(f"Percent of host genome found: {outside_args/2011548}\n")
        
        #output.write(f"Area of args that were inserted: {area_of_args}\n")
        #output.write(f"bp area that filter found within ARGs: {within_args}\n")
        #output.write(f"total bp area of the args that we did find: {all_area_of_found_args}\n")
        #output.write(f"bp area that filter found outside ARGs: {outside_args}\n")        

        
        output.write(f"Of args found, the percent of area found within them: {within_args / all_area_of_found_args}\n")
        output.write(f"Of args found, percent that has >= 90% of area found: {num_over_90 / len(found_args)}")
        
        #output.write("Found args, and the bp of them found::\n")
        #for arg in found_args:
        #    output.write(f"{arg},{found_args[arg]}\n")
