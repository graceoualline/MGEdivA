# if a file is over 200 kbp, this will sqpit the query into 200kbp files
# right now, even if it is 200,001 bp it will get split, will worry about this later

from Bio import SeqIO
from pathlib import Path

def split_query(input_fasta, chunk_size=200_000):
    out_dir = input_fasta.split(".")[0] + "_chunks"
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    records = list(SeqIO.parse(input_fasta, "fasta"))
    chunk_files = []

    for record in records:
        seq = record.seq
        for i in range(0, len(seq), chunk_size):
            chunk_seq = seq[i:i+chunk_size]
            chunk_id = f"{record.id}_chunk_{i}"
            chunk_file = Path(out_dir) / f"{chunk_id}.fasta"
            with open(chunk_file, "w") as out_f:
                out_f.write(f">{chunk_id}\n{chunk_seq}\n")
            chunk_files.append(str(chunk_file))

    return chunk_files