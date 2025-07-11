import csv
import sys
import pandas as pd

def get_coor(AMRFinder):
    coor = []
    with open(AMRFinder, 'r') as infile:
        # Process each line in the file
        for line in infile:
            columns = line.strip().split('\t')
            start, end = columns[2], columns[3]
            if not start.isnumeric():
                continue
            coor.append((start, end))
    return coor

def region_overlap(r1, r2):
    return max(0, min(float(r1[1]), float(r2[1])) - max(float(r1[0]), float(r2[0])))

def analyze(true_args, coor):
    found_args = {}
    total_overlap_area = 0
    total_found_area = 0

    for i, c in enumerate(coor): 
        #print(i)
        total_found_area += int(c[1])-int(c[0])
        for j in range(len(true_args)):
            start = true_args.iloc[j, 1]
            end = true_args.iloc[j, 2]
            overlap = region_overlap(c, (start, end))
            if overlap > 0:
                arg = true_args.iloc[j, 0]
                if arg in found_args: found_args[arg] += overlap
                else: found_args[arg] = overlap
                #this is assuming AMR doesnt return stuff that overlaps
                total_overlap_area += overlap
    return found_args, total_overlap_area, total_found_area - total_overlap_area

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python3 analyze_AMRFinder.py <AMR_finder_output> <.csv file of locations and names> <output file of stats>")
        sys.exit(1)
    # Example usage
    AMRFinder = sys.argv[1]
    locations = sys.argv[2]
    output = sys.argv[3]

    locations = pd.read_csv(locations)
    #locations = locations[locations.iloc[:, 1] <= 11500000]  # REMOVE THIS IF WE END UP DOING MORE THAN 114 CHUNKS, filtered for args within the first 114 chunks
    print(locations[:3])

    coordinates = get_coor(AMRFinder)
    print(f"{len(coordinates)} coordinates extracted")
    #print(coordinates[:3])

    # Write the overlapping plasmids to a .txt file
    found_args, within_args, outside_args = analyze(locations, coordinates)

    arg_regions = locations
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
'''
    with open(output_file, 'w') as output:
        output.write(f"Number of ARGs found: {len(found)}\n")
        output.write(f"bp area found within ARGs: {area_correct}\n")
        output.write(f"bp area found outside ARGs: {area_incorrect}\n")
        output.write(f"percentage of args found: {len(found) / 10000}\n")
        output.write(f"percentage of arg area found: {area_correct / 10230387}\n")
        output.write(f"percentage of non-arg found : {area_incorrect / 1269613}\n")
        output.write(f"Found args:\n")
        for arg in found:
            output.write(f"{arg}\n")
    
    print("Change was saved successfully")
  '''  