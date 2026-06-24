import os, glob, polars as pl, numpy as np
from pathlib import Path

def parse_phages(fasta_files):
    print(f"[Phase 2] Parsing {len(fasta_files)} FASTA sequence file(s)...")
    records = []
    for p in fasta_files:
        v_id = Path(p).stem
        g, c, t = 0, 0, 0
        with open(p, 'rb') as f:
            for l in f:
                if l.startswith(b">"): continue
                cl = l.strip().upper(); t += len(cl); g += cl.count(b"G"); c += cl.count(b"C")
        records.append({"Viral_Identity": v_id, "Genome_Length_BP": t, "True_GC_Content_Pct": round(((g+c)/t*100), 1) if t > 0 else 0.0})
    return pl.DataFrame(records)

def scan_multiple_motifs(fasta_files, motifs=["CGATC", "GATC", "CCGG"]):
    print(f"[Phase 3] Scanning sequence matrices for target motifs: {motifs}...")
    records = []
    for p in fasta_files:
        v_id = Path(p).stem
        with open(p, 'rb') as f: 
            r_seq = b"".join([l.strip().upper() for l in f if not l.startswith(b">")])
        
        entry = {"Viral_Identity": v_id}
        for m in motifs:
            entry[f"Motif_{m}_Count"] = r_seq.count(m.encode('utf-8'))
        records.append(entry)
    return pl.DataFrame(records)

def train_model(m_df):
    print("\n--- Running Phase 4: Predictive Analytics Engine ---")
    numeric_cols = [c for c in m_df.columns if c != "Viral_Identity"]
    
    unique_rows = m_df.select(numeric_cols).unique()
    if unique_rows.height < 3:
        print("⚠️ Insufficient biological variance or too few unique strains to map genomic correlations.")
        return

    print("✅ Biological variance detected! Calculating genomic correlation matrix...")
    corr_matrix = m_df.select(numeric_cols).corr()
    print(corr_matrix)

def run_genomics_pipeline():
    print("🚀 Initializing Cross-Disciplinary Genomics Framework...")
    
    # 🎯 Target the exact folder where the fetcher saved the 10 files
    data_dir = Path("ncbi_dataset/data")
    f_files = [str(p) for p in data_dir.rglob("*") if p.suffix in [".fna", ".fasta"]]

    if len(f_files) >= 2:
        r_df = parse_phages(f_files)
        w_df = scan_multiple_motifs(f_files)
        f_matrix = r_df.join(w_df, on="Viral_Identity", how="inner")
        print(f"\n--- Processed Genetic Report ---\n{r_df}\n\n--- Phage Weapon Discovery Report ---\n{w_df}")
        train_model(f_matrix)
        f_matrix.write_csv("klebsiella_phage_matrix_summary.csv")
        print("\nScan complete! Results successfully saved to 'klebsiella_phage_matrix_summary.csv'.")
    else:
        print("🧬 No local raw datasets detected in ncbi_dataset/data. Running pipeline with simulated dataset...")
        r_df = pl.DataFrame([
            {"Viral_Identity":"Phage_Sweeny", "Genome_Length_BP": 41119, "True_GC_Content_Pct": 52.4},
            {"Viral_Identity":"Phage_KpGranit", "Genome_Length_BP": 45200, "True_GC_Content_Pct": 49.1}
        ])
        w_df = pl.DataFrame([
            {"Viral_Identity":"Phage_Sweeny", "Motif_CGATC_Count": 1, "Motif_GATC_Count": 12, "Motif_CCGG_Count": 4},
            {"Viral_Identity":"Phage_KpGranit", "Motif_CGATC_Count": 4, "Motif_GATC_Count": 22, "Motif_CCGG_Count": 9}
        ])
        f_matrix = r_df.join(w_df, on="Viral_Identity", how="inner")
        print(f_matrix)

if __name__ == "__main__": 
    run_genomics_pipeline()
