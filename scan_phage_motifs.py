import polars as pl

print("====================================================")
print("  KLEBSIELLA PHAGE DEPOLYMERASE MOTIF SCANNER v2.0  ")
print("====================================================\n")

phage_names = []
sequences = []

with open("superbug_phages.fasta", "r") as file:
    current_seq = []
    for line in file:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq))
                current_seq = []
            phage_names.append(line[1:])
        else:
            current_seq.append(line)
    if current_seq:
        sequences.append("".join(current_seq))

# Multi-motif profiles for complex structural scanning
motifs = {
    "Tail_Hinge_GGTA": "GGTA",
    "Capsule_Binder_GACG": "GACG",
    "Anchor_CGATC": "CGATC"
}

results = []
for name, seq in zip(phage_names, sequences):
    seq = seq.upper()
    ggta_hits = seq.count(motifs["Tail_Hinge_GGTA"])
    gacg_hits = seq.count(motifs["Capsule_Binder_GACG"])
    cgatc_hits = seq.count(motifs["Anchor_CGATC"])
    weapon_score = ggta_hits + gacg_hits + cgatc_hits
    
    results.append({
        "Phage_Identity": name[:30] + "...", # Shorten name for clean scannable layout
        "Hinge_Matches": ggta_hits,
        "Binder_Matches": gacg_hits,
        "Anchor_Matches": cgatc_hits,
        "Total_Weapon_Score": weapon_score
    })

analysis_df = pl.DataFrame(results)
ranked_phages = analysis_df.sort("Total_Weapon_Score", descending=True)

print("--- Real-World Phage Weapon Analysis Complete ---")
print(ranked_phages)
print("-------------------------------------------------")
