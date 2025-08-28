# this code will intake a gtdb multifasta file or fasta file of your choice
# and give you a fasta database that you can use for blat
from Bio import SeqIO
import sys
import os
import subprocess
import time

def write_group(group_seqs, index):
        # if somehow the list is empty, skip
        if not group_seqs:
            return
        out_path = os.path.join(output_dir, f"split_{index}.fasta")
        with open(out_path, "w") as handle:
            SeqIO.write(group_seqs, handle, "fasta")
        print(f"Wrote {out_path}")

def progress(processed_seqs, total_seqs, start_time, file_index):
    # Progress indicator
    if processed_seqs % 10000 == 0 or processed_seqs == total_seqs:
        elapsed = time.time() - start_time
        percent = (processed_seqs / total_seqs) * 100
        rate = processed_seqs / elapsed if elapsed > 0 else 0
        eta = (total_seqs - processed_seqs) / rate if rate > 0 else 0
        
        print(f"Progress: {processed_seqs:,}/{total_seqs:,} sequences ({percent:.1f}%) | "
                f"Rate: {rate:.0f} seq/s | ETA: {eta:.0f}s | Files created: {file_index-1}")

def count_sequences(input_fasta):
    print(f"Counting sequences in {input_fasta}...")
    count = 0
    for _ in SeqIO.parse(input_fasta, "fasta"):
        count += 1
        if count % 10000 == 0: print("Current count:", count)
    return count

def count_sequences_in_directory(directory):
    """Count total sequences in all FASTA files in a directory"""
    if not os.path.exists(directory):
        return 0
    
    total_count = 0
    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
    
    for fasta_file in fasta_files:
        fa = os.path.join(directory, fasta_file)
        total_count += count_sequences(fa)
    return total_count

def check_step1_completion(input_fasta, output_dir, original_count):
    """Check if step 1 is complete by comparing sequence counts"""
    if not os.path.exists(output_dir):
        return False
    
    # Count sequences in split files
    split_count = count_sequences_in_directory(output_dir)
    
    is_complete = original_count == split_count and split_count > 0
    
    print(f"Step 1 status check:")
    print(f"  Original file sequences: {original_count:,}")
    print(f"  Split files sequences: {split_count:,}")
    print(f"  Step 1 complete: {'Yes' if is_complete else 'No'}")
    
    return is_complete

def split_fasta_by_bp(input_fasta, output_dir, max_bp):
    #make dir if it doesnt exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Count total sequences for progress tracking
    total_seqs = count_sequences(input_fasta)

    complete = check_step1_completion(input_fasta, output_dir, total_seqs)
    if complete: return 1

    group = []
    group_bp = 0
    file_index = 1
    processed_seqs = 0
    start_time = time.time()

    print("\n" + "="*60)
    print("STEP 1: Splitting FASTA file by base pairs")
    print("="*60)

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_len = len(record.seq)

        processed_seqs += 1
        progress(processed_seqs, total_seqs, start_time, file_index)
        
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

    total_time = time.time() - start_time
    print(f"\nStep 1 completed in {total_time:.1f}s")
    print(f"Created {file_index} FASTA files")

def convert_2bit(output_name):
    input_directory = output_name + "_fa"
    output_directory = output_name + "_2bit"
    if not os.path.exists(output_directory):
            os.makedirs(output_directory)

    # Get list of files to process
    fasta_files = [f for f in os.listdir(input_directory) 
                   if os.path.isfile(os.path.join(input_directory, f))]
    total_files = len(fasta_files)

    print("\n" + "="*60)
    print("STEP 2: Converting FASTA files to 2bit format")
    print("="*60)
    print(f"Total files to convert: {total_files}")

    file_type = "2bit"
    start_time = time.time()

    for i, input_file in enumerate(fasta_files, 1):
        file_name, file_extension = os.path.splitext(input_file)
        
        # Call your Python script with the input and output file names
        input_path = os.path.join(input_directory, input_file)
        output_path = os.path.join(output_directory, f"{file_name}.{file_type}")
        # if the output file already exists, we don't need to make it again
        if os.path.exists(output_path): continue
        # Progress indicator
        percent = (i / total_files) * 100
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_files - i) / rate if rate > 0 else 0
        
        print(f"Progress: {i}/{total_files} files ({percent:.1f}%) | "
              f"Rate: {rate:.1f} files/s | ETA: {eta:.0f}s")

        try:
            result = subprocess.run(["faToTwoBit", input_path, output_path], 
                                  capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR converting {input_file}: {e}")
            print(f"stderr: {e.stderr}")
        except FileNotFoundError:
            print("ERROR: faToTwoBit command not found. Please ensure it's installed and in your PATH.")
            break

    total_time = time.time() - start_time
    print(f"\nStep 2 completed in {total_time:.1f}s")
    print(f"Converted {total_files} files to 2bit format")

def make_ooc(output_name):
    input_directory = output_name + "_2bit"
    bit_files = [f for f in os.listdir(input_directory) 
                 if f.endswith('.2bit') and os.path.isfile(os.path.join(input_directory, f))]

    total_files = len(bit_files)
    
    print("\n" + "="*60)
    print("STEP 3: Generating OOC files from 2bit files")
    print("="*60)
    print(f"Total 2bit files to process: {total_files}")

    if total_files == 0:
        print("ERROR: No 2bit files found in directory:", input_directory)
        return

    start_time = time.time()

    for i, bit_file in enumerate(bit_files, 1):
        # Generate output filename (replace .2bit with .ooc)
        base_name = bit_file.replace('.2bit', '')
        # Check if the file is a regular file
        input_path = os.path.join(input_directory, bit_file)
        output_file = f"{base_name}.ooc"
        output_path = os.path.join(input_directory, output_file)
        # if the output file already exists, we don't need to make it again
        if os.path.exists(output_path): continue

        percent = (i / total_files) * 100
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total_files - i) / rate if rate > 0 else 0
        print(f"Progress: {i}/{total_files} files ({percent:.1f}%) | "
              f"Rate: {rate:.1f} files/s | ETA: {eta:.0f}s")

        try:
            cmd = ["blat", input_path, "/dev/null", "/dev/null", f"-makeOoc={output_path}", "-repMatch=1024"]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"ERROR generating OOC for {bit_file}: {e}")
            print(f"stderr: {e.stderr}")
        except FileNotFoundError:
            print("ERROR: blat command not found. Please ensure it's installed and in your PATH.")
            break
    total_time = time.time() - start_time
    print(f"\nStep 3 completed in {total_time:.1f}s")
    print(f"Generated OOC files for {total_files} 2bit files")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python make_blat_db.py <input.fasta> <output_dir> <max_bil_bp>")
        sys.exit(1)
    
    input_fasta = sys.argv[1] #full_database = "/usr1/shared/gtdb_combined.fa"
    output_dir = sys.argv[2] + "_fa"
    max_bp = int(sys.argv[3]) * 1000000000

    split_fasta_by_bp(input_fasta, output_dir, max_bp)

    print("\nFinished splitting the query file. Converting to 2bit.")

    convert_2bit(sys.argv[2])

    print("\nFinished 2bit conversion. Generating OOC files.")

    make_ooc(sys.argv[2])

    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)
    

# python3 /usr1/gouallin/blat/blat_pipeline/make_blat_db.py /usr1/shared/gtdb_combined.fa gtdb_smart_split_fa 2
