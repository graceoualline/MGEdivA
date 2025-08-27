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

def split_fasta_by_bp(input_fasta, output_dir, max_bp):
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
            write_group(group, file_index)
            file_index += 1
            group = [record]
            group_bp = seq_len
        else:
            group.append(record)
            group_bp += seq_len

    # Write the final group
    if group:
        write_group(group, file_index)

def convert_2bit(output_name):
    input_directory = output_name + "fa"
    output_directory = output_name + "_2bit"
    if not os.path.exists(output_directory):
            os.makedirs(output_directory)
    file_type = "2bit"
    for input_file in os.listdir(input_directory):
        # Check if the file is a regular file
        if os.path.isfile(os.path.join(input_directory, input_file)):
            # Extract the file name without the extension
            file_name, file_extension = os.path.splitext(input_file)
            
            # Call your Python script with the input and output file names
            input_path = os.path.join(input_directory, input_file)
            output_path = os.path.join(output_directory, f"{file_name}_output.{file_type}")
            print("Running", input_path, output_path)
            subprocess.run(["faToTwoBit", input_path, output_path])

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python make_blat_db.py <input.fasta> <output_dir> <max_bil_bp>")
        sys.exit(1)
    
    input_fasta = sys.argv[1] #full_database = "/usr1/shared/gtdb_combined.fa"
    output_dir = sys.argv[2] + "_fa"
    max_bp = int(sys.argv[4]) * 1000000000

    split_fasta_by_bp(input_fasta, output_dir, max_bp)

    print("Finished splitting the query file. Converting to 2bit.")

    convert_2bit(sys.argv[2])
    

# python3 /usr1/gouallin/blat/blat_pipeline/split_fasta_by_bp.py /usr1/shared/gtdb_combined.fa gtdb_smart_split_fa 5 6
