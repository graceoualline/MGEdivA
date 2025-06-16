#I willt ake all of the blat 2bit files and make ooc files for them

import os
import sys
import subprocess

input_files = [f"gtdb1k_part_{i}_output.2bit" for i in range(1, 1001)]
output_names = [f"part{i}_1k_11.ooc" for i in range(1, 1001)]

for input_file, output_file in zip(input_files, output_names):
    # Check if the file is a regular file
    print("running:", "blat", input_file, "/dev/null", "/dev/null", f"-makeOoc={output_file}", "-repMatch=1024")
    subprocess.run(["blat", input_file, "/dev/null", "/dev/null", f"-makeOoc={output_file}", "-repMatch=1024"])