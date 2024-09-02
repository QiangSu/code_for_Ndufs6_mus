import pandas as pd
import os

# Define the file names
files = [
    "./ADULT1/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_ADULT1.csv",
    "./ADULT2/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_ADULT2.csv",
    "./ADULT3/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_ADULT3.csv",
    "./KO1/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KO1.csv",
    "./KO2/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KO2.csv",
    "./KO3/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KO3.csv",
    "./KO_AAV1/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KOAAV1.csv",
    "./KO_AAV2/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KOAAV2.csv",
    "./KO_AAV3/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_KOAAV3.csv",
    "./NB1/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_NB1.csv",
    "./NB2/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_NB2.csv",
    "./NB3/sorted_genes.fpkm_tracking_cufflinks_result_nonbiascorrect_NB3.csv",

]

# Initialize an empty DataFrame for the combined data
combined_df = pd.DataFrame()

for file in files:
    # Read the CSV file with tab separation
    df = pd.read_csv(file, delimiter='\t')

    # Initialize the combined DataFrame with required columns only if it is empty
    if combined_df.empty:
        combined_df = df[['tracking_id', 'gene_id', 'gene_short_name', 'length']].copy()

    # Extract the sample name from the file name
    base_name = os.path.basename(file)
    sample_name_parts = base_name.split('_')[6:7]
    sample_name = '_'.join(sample_name_parts).replace('.csv', '')

    # Append the coverage column to the combined DataFrame
    combined_df[sample_name + '_coverage'] = df['coverage']

# Save the combined DataFrame to a new CSV file
combined_df.to_csv('mus_cufflinks_nobiascorrect_gene_combined_coverage.csv', index=False)
