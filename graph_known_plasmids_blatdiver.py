import matplotlib.pyplot as plt
import matplotlib.patches as patches
from extract_blatdiver_coor import *
import sys
from matplotlib.patches import Patch


def extract_special_blatdiver_coor(input):
    query_coor = []
    with open(input, "r") as f:
        for line in f:
            #print(line)
            cols = line.strip().split("\t")
            #print("cols", cols)
            #print("col 0", cols)
            if len(cols) == 0 or cols[0] == "Q name":  # Skip header if it exists
                continue
            #print(cols[11])
            q_start = int(cols[2])  # Q start (0-based) 
            q_end = int(cols[3])    # Q end
            # Write coordinates to the output file
            query_coor.append((q_start,q_end))

    print(f"Coordinates extracted from {input}")
    return query_coor

def plot_graph(coor, gene_length, title, output_file, plasmid_loc):
    skani_color =  "red"
    plasmid_color = "yellow"
    edge_color = "purple"
    #genomad_color =  (1/255, 207/255, 13/255)
    #mefinder_color =  (1/255, 40/255, 218/255)

    # Create a figure and axis
    fig, axes = plt.subplots(5, 1, figsize=(9, 2), sharex=True, gridspec_kw={'hspace': 2}) #fig, ax = plt.subplots(figsize=(9, 0.5))

    #print(coor)
    plasmid_data = [coor[i] for i in range(len(coor)) if (coor[i][0] >= plasmid_loc[0] and coor[i][1] <= plasmid_loc[1])]
    print("plasmid data num", len(plasmid_data))
    skani_data = [coor[i] for i in range(len(coor)) if coor[i][1] < plasmid_loc[0] or coor[i][0] > plasmid_loc[1]]
    print("skani data num", len(skani_data))
    edge_data = list(set(coor) - (set(plasmid_data).union(set(skani_data))))
    print("edge data num", len(edge_data))
    known_location = [plasmid_loc]

    subplots_info = [
        (title, coor, [(edge_color, edge_data), (skani_color, skani_data), (plasmid_color, plasmid_data)]),
        ("Blatdiver Hits Outside Plasmid", skani_data, [(skani_color, skani_data)]),
        ("In Plasmid", plasmid_data, [(plasmid_color, plasmid_data)]),
        ("Where Plasmid is", known_location, [("green", known_location)]),
        ("Overlapping In and Out of Plasmid", edge_data, [(edge_color, edge_data)])
        
    ]
    for ax, (subplot_title, data, color_data) in zip(axes, subplots_info):
        for color, data_points in color_data:
            for start, end in data_points:
                rect = patches.Rectangle((start, 0), end - start, 1, linewidth=1, edgecolor=color, facecolor=color, alpha=0.5)
                ax.add_patch(rect)

        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_xlim(0, int(gene_length))
        ax.set_title(subplot_title, fontsize=10, pad=10)

    # Common X-axis settings
    axes[-1].set_xticks(range(0, gene_length, gene_length // 10))
    axes[-1].set_xlabel("Genomic Position")

    # Save the figure
    plt.savefig(output_file, bbox_inches="tight", dpi=300)
    print("Figure saved to", output_file)
'''
    for start3, end3 in edge_data:
        rect3 = patches.Rectangle((start3, 0), end3 - start3, 1, linewidth=1, edgecolor= edge_color, facecolor=edge_color, alpha=0.1)
        ax.add_patch(rect3)
    for start1, end1 in skani_data:
        rect1 = patches.Rectangle((start1, 0), end1 - start1, 1, linewidth=1, edgecolor= skani_color, facecolor=skani_color, alpha=0.5)
        ax.add_patch(rect1)
    for start2, end2 in plasmid_data:
        rect2 = patches.Rectangle((start2, 0), end2 - start2, 1, linewidth=1, edgecolor= plasmid_color, facecolor=plasmid_color, alpha=0.5)
        ax.add_patch(rect2)

    # Hide y-axis
    ax.set_yticks([])
    ax.set_yticklabels([])
    maximum = int(gene_length)


    # Add labels
    ax.set_xticks(range(0,  maximum, maximum//10))
    ax.set_xlabel('Genomic Position')


    plt.title(title)

    


    #plt.show()
    # or save the plot
    plt.savefig(output_file, bbox_inches='tight', dpi=300)
    print("Figure saved to", output_file)
'''

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python3 graph_known_plasmids_blatdiver.py <blatdiver output> <query gene length> <figure title> <png output file name> <is filtered (T/F)>")
        sys.exit(1)

    # Usage example:
    blatdiver = sys.argv[1]
    gene_length = int(sys.argv[2])
    title = sys.argv[3]
    output_file = sys.argv[4]
    filtered = sys.argv[5]
    if filtered == "True":
        coor = extract_special_blatdiver_coor(blatdiver) #this is for overlap filtered files
    else:
        coor = extract_blatdiver_coor(blatdiver) #this is for non-overlap filtered files
    
    plasmid_loc = blatdiver.split(",")
    #print("plasmid_loc", plasmid_loc)
    #plasmid_loc = (10507,123111)
    plasmid_loc = (int(plasmid_loc[1]), int(plasmid_loc[2]))
    plot_graph(coor, gene_length, title, output_file, plasmid_loc)

'''
input = "NC_016607.1,2488,25624,CP073938.1_filtered.tsv"
coor = extract_special_blatdiver_coor(input)
plasmid_loc = (2488,25624)
gene_length = 26424
title = "NC_016607.1,2488,25624,CP073938.1_filtered"
output_file = "NC_016607.1,2488,25624,CP073938.1_filtered.png"
plot_graph(coor, gene_length, title, output_file, plasmid_loc)

'''