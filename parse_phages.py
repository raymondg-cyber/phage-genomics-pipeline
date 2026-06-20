import polars as pl

print("Opening FASTA file and reading raw viral sequence strings...")

phage_names = []
sequences = []

# Open and parse the plain text FASTA file line-by-line (perfect for 8GB Macs!)
with open("superbug_phages.fasta", "r") as file:
    current_seq = []
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            # If we hit a new phage header, save the previous one first
            if current_seq:
                sequences.append("".join(current_seq))
                current_seq = []
            phage_names.append(line[1:])  # Strip out the '>' character
        else:
            current_seq.append(line)
    # Don't forget to append the final sequence in the file
    if current_seq:
        sequences.append("".join(current_seq))

# Calculate raw text stats for each virus
gc_percentages = []
sequence_lengths = []

for seq in sequences:
    # Convert all characters to uppercase to prevent counting bugs
    seq = seq.upper()
    seq_length = len(seq)
    
    # Count the specific structural targeting letters (G and C)
    g_count = seq.count("G")
    c_count = seq.count("C")
    
    # Calculate the definitive GC percentage
    gc_pct = ((g_count + c_count) / seq_length) * 100
    
    gc_percentages.append(round(gc_pct, 2))
    sequence_lengths.append(seq_length)

# Port our raw text extraction into a clean Polars Dataframe
results_df = pl.DataFrame({
    "Viral_Identity": phage_names,
    "Genome_Length_BP": pl.Series(sequence_lengths, dtype=pl.Int32),
    "True_GC_Content_Pct": pl.Series(gc_percentages, dtype=pl.Float32)
})

# Save our molecular extraction report
results_df.write_csv("phage_genetic_report.tsv", separator="\t")

print("\n--- Processed Genetic Report ---")
print(results_df)
print("--------------------------------")
print("Success! Viral sequences successfully parsed and quantified.")
