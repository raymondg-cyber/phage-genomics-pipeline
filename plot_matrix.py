import os
import polars as pl
import matplotlib.pyplot as plt
import seaborn as sns

def generate_heatmap(csv_path="klebsiella_phage_matrix_summary.csv", output_path="correlation_heatmap.png"):
    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found. Run the pipeline first!")
        return

    # Load data using Polars and compute the Pearson correlation matrix
    df = pl.read_csv(csv_path)
    numeric_df = df.select([
        "Genome_Length_BP", "True_GC_Content_Pct", 
        "Motif_CGATC_Count", "Motif_GATC_Count", "Motif_CCGG_Count"
    ])
    corr_matrix = numeric_df.corr().to_pandas()
    corr_matrix.index = numeric_df.columns

    # Set up professional publication-style plotting aesthetics
    plt.figure(figsize=(8, 6), dpi=300)
    sns.set_theme(style="white")

    # Generate custom diverging colormap frequently used in genomics
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Plot the matrix
    sns.heatmap(
        corr_matrix, 
        annot=True, 
        fmt=".3f", 
        cmap=cmap, 
        vmax=1.0, 
        vmin=0.5, 
        square=True, 
        linewidths=.5, 
        cbar_kws={"shrink": .75, "label": "Pearson Correlation Coefficient (r)"}
    )

    plt.title("Genomic Correlation Matrix: Klebsiella Phage vs. Host Attributes", fontsize=12, fontweight="bold", pad=15)
    plt.tight_layout()

    # Save image asset for GitHub README rendering
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Success: Publication-quality heatmap saved to '{output_path}'")

