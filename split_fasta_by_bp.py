# this code will intake a gtdb multifasta file or fasta file of your choice
# and give you a fasta database that you can use for blat
from Bio import SeqIO
import sys
import os

def write_group(group_seqs, index):
        # if somehow the list is empty, skip
        if not group_seqs:
            return
        out_path = os.path.join(output_dir, f"split_{index}.fasta")
        with open(out_path, "w") as handle:
            SeqIO.write(group_seqs, handle, "fasta")
        print(f"Wrote {out_path}")

def split_fasta_by_bp(input_fasta, output_dir, min_bp, max_bp):
    #make dir if it doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    group = []
    group_bp = 0
    file_index = 1

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)
        

        # If one sequence is longer than max_bp, write it by itself
        if seq_len > max_bp:
            print(f"Warning: {record.id} is {seq_len} bp, longer than max ({max_bp} bp). Saving in its own file.")
            write_group([record], file_index)
            file_index += 1
            continue

        # If adding this record would exceed max_bp, flush and start new group
        if group_bp + seq_len > max_bp:
            #aka if min_bp < group_bp < max_bp
            if group_bp >= min_bp:
                write_group(group, file_index)
                file_index += 1
                group = [record]
                group_bp = seq_len
            else:
                # Try to squeeze in more until we reach min_bp
                group.append(record)
                group_bp += seq_len
        else:
            group.append(record)
            group_bp += seq_len

    # Write the final group
    if group:
        write_group(group, file_index)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python split_fasta_by_bp.py <input.fasta> <output_dir> <min_gb> <max_gb>")
        sys.exit(1)
    # To figure out how big you want each file to be,
    # calculate based on the maximum number of threads you
    # plan to use at any time, and divide the total avaliable RAM you
    # have by the thread amount. For example
    # I use 46 threads at most when running blatdiver
    # I have about 350 GB of avaliable RAM
    # 380 GB / 46 threads = ~8.2 GB of RAM per thread
    # So I will want my file side to be maxmimum 6, minimum 5 GB large
    
    input_fasta = sys.argv[1] #full_database = "/usr1/shared/gtdb_combined.fa"
    output_dir = sys.argv[2]
    # 4e9 bp = 1 GB
    min_bp = int(sys.argv[3]) * 1000000000
    max_bp = int(sys.argv[4]) * 1000000000

    split_fasta_by_bp(input_fasta, output_dir, min_bp, max_bp)

# python3 /usr1/gouallin/blat/blat_pipeline/split_fasta_by_bp.py /usr1/shared/gtdb_combined.fa gtdb_smart_split_fa 5 6
