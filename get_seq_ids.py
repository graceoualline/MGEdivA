#will extract all of the sequence ids from teh 2bit files and put them into .txt files
import os
import subprocess

# Generate file names
input_files = [f"gtdb1k_part_{i}_output.2bit" for i in range(1, 1001)]
output_files = [f"part_{i}_info.txt" for i in range(1, 1001)]
final_output = "sequence_id_to_2bitfile_mapping.txt"  # Consolidated output file

with open(final_output, "w") as fout:
    for inputf, output in zip(input_files, output_files):
        two_bit_path = f"/usr1/shared/gtdb_split_2bit_1k/{inputf}"
        info_path = f"2bit_info/{output}"

        #print(f"Processing {two_bit_path} → {info_path}")

        # Run twoBitInfo to get sequence info
        print("running")
        print("twoBitInfo", two_bit_path, info_path)
        subprocess.run(["twoBitInfo", two_bit_path, info_path])
        

        # Read the sequence IDs and write with file name
        with open(info_path, "r") as f:
            for line in f:
                seq_id = line.split()[0]
                fout.write(f"{seq_id}\t{inputf}\n")
