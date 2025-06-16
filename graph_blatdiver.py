import matplotlib.pyplot as plt
import matplotlib.patches as patches
from extract_blatdiver_coor import *
import sys

def plot_graph(coor, gene_length, title, output_file, color = "red"):
    skani_color =  color
    #genomad_color =  (1/255, 207/255, 13/255)
    #mefinder_color =  (1/255, 40/255, 218/255)

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(9, 0.5))

    max_s = 0
    max_m = 0
    max_g = 0


    skani_data = coor

    for start1, end1 in skani_data:
        rect1 = patches.Rectangle((start1, 0), end1 - start1, 1, linewidth=1, edgecolor= skani_color, facecolor=skani_color, alpha=0.5)
        ax.add_patch(rect1)

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

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 5:
        print("Usage: python3 graph_blatdiver.py <blatdiver output> <query gene length> <figure title> <png output file name>")
        sys.exit(1)

    # Usage example:
    blatdiver = sys.argv[1]
    gene_length = sys.argv[2]
    title = sys.argv[3]
    output_file = sys.argv[4]
    coor = extract_blatdiver_coor(blatdiver)
    plot_graph(coor, gene_length, title, output_file)
