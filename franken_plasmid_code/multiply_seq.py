#I need a certain length for the host genome to get the result I want
# was used to multiplle the host genome 66 times so it is long enough for the franken plasmid
#ChatGPT
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# === CONFIGURATION ===
input_fasta = "CP124955.1_host_genome.fasta"      # <- change this
output_fasta = "CP124955.1_host_genome_multx66.fasta"
repeat_count = 66

# === READ THE SEQUENCE ===
record = next(SeqIO.parse(input_fasta, "fasta"))
print(f"Read sequence: {record.id} (length: {len(record.seq)})")

# === MULTIPLY IT ===
multiplied_seq = record.seq * repeat_count
print(f"Multiplied length: {len(multiplied_seq)}")

# === CREATE A NEW RECORD ===
new_record = SeqRecord(
    multiplied_seq,
    id=f"{record.id}_x{repeat_count}",
    description=f"{record.description} multiplied {repeat_count} times"
)

# === WRITE TO FILE ===
SeqIO.write(new_record, output_fasta, "fasta")
print(f"Written multiplied sequence to {output_fasta}")