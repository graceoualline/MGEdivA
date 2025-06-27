#will select random args that we will use for franken args
#Chat GPT
import random
from Bio import SeqIO

# === CONFIGURATION ===
input_fasta = "/usr1/gouallin/mobile_ctrl/args/nucleotide_fasta_protein_homolog_model_variants.fasta"      # <- change this
output_fasta = "rand_1000_args.fasta"
num_sequences = 1000

# === LOAD ALL SEQUENCES ===
print("Reading sequences from file...")
records = list(SeqIO.parse(input_fasta, "fasta"))
total = len(records)

if total < num_sequences:
    raise ValueError(f"File only contains {total} sequences, but you requested {num_sequences}.")

# === RANDOMLY SAMPLE ===
print(f"Randomly selecting {num_sequences} sequences from {total} total...")
subset = random.sample(records, num_sequences)

# === WRITE TO FILE ===
print(f"Writing to {output_fasta}...")
SeqIO.write(subset, output_fasta, "fasta")

print("Done!")
