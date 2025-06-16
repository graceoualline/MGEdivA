#ai generated


def extract_blatdiver_coor(input):
    query_coor = []
    with open(input, "r") as f:
        for line in f:
            #print(line)
            cols = line.strip().split("\t")
            #print("cols", cols)
            #print("col 0", cols)
            if len(cols) == 0 or cols[0] == "Q name":  # Skip header if it exists
                continue
            if not (cols[0]).isnumeric():  # Ensure valid format
                continue
            #print(cols[11])
            q_start = int(cols[11])  # Q start (0-based) 
            q_end = int(cols[12])    # Q end
            # Write coordinates to the output file
            query_coor.append((q_start,q_end))

    print(f"Coordinates extracted from {input}")
    return query_coor
