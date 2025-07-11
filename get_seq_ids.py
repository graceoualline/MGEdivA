#will extract all of the sequence ids from teh 2bit files and put them into .txt files
import os
import subprocess
import tempfile

# Generate file names
def main():
    input_files = [f"split_{i}_output.2bit" for i in range(1, 138)]
    final_output = "sequence_id_to_2bit_mapping_smart_split.txt"  # Consolidated output file

    with open(final_output, "w") as fout:
        for inputf in input_files:
            two_bit_path = f"/usr1/gouallin/blat/gtdb_smart_split_2bit/{inputf}"

            #this will make a temp file
            with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
                temp_path = temp_file.name

            # Run twoBitInfo to get sequence info
            try:
                print("running")
                print("twoBitInfo", two_bit_path, temp_path)
                result = subprocess.run(["twoBitInfo", two_bit_path, temp_path],
                                        capture_output = True,
                                        text = True,
                                        check = True)
            
                # Read the sequence IDs and write with file name
                with open(temp_path, "r") as f:
                    for line in f:
                        seq_id = line.split()[0]
                        fout.write(f"{seq_id}\t{inputf}\n")
            finally:
                os.remove(temp_path)
                #assert(False)
