import pysam
import gzip
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import os
import argparse
import seaborn as sns
import sys

print("Sample Name:", sys.argv[1])

sample_folder = f'/home/disk/chengrui/workspace/spatial_transcriptomics/data/summarized_data/{sys.argv[1]}'
os.chdir(sample_folder)

# Initialize dictionaries to store gene, barcode, and UMI data
new_gene, old_gene = set(), set()
new_gene_umi, old_gene_umi = {}, {}
bc_new_gene, bc_old_gene = {}, {}
bc_new_gene_umi, bc_old_gene_umi = {}, {}

# Load cell barcodes from the filtered barcode file
with gzip.open(f'{sample_folder}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz', 'rt') as f:
    cell_barcode = set(line.strip() for line in f)

def process_gene_data(barcode, geneid, umi, gene, gene_dict, gene_umi_dict, bc_gene_dict, bc_gene_umi_dict):
    """
    Process the gene data for new or old genes, updating dictionaries for gene, barcode, and UMI.
    """
    gene_name = gene.rstrip('-C') if gene.endswith('--C') else gene.rstrip('--T')
    if gene.endswith('--T') or gene.endswith('--C'):
        gene_dict.add(gene_name)
        gene_umi_dict.setdefault(gene_name, set()).add(f'{barcode}#{umi}')
        bc_gene_dict.setdefault(barcode, set()).add(gene_name)
        bc_gene_umi_dict.setdefault(barcode, set()).add(f'{gene_name}#{umi}')

# Read BAM file and extract barcode, gene, and UMI information
with pysam.AlignmentFile(f'{sample_folder}/Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam', 'rb') as f:
    for r in f.fetch():
        barcode, geneid, umi = r.get_tag("CB").replace("_",""), r.get_tag("GX"), r.get_tag("UB")
        
        # Process only valid reads
        if barcode != "-" and geneid != "-" and umi != "-" and barcode in cell_barcode:
            try:
                gene = r.get_tag('GN')
                # Process genes
                if gene.endswith('--T'):
                    process_gene_data(barcode, geneid, umi, gene, old_gene, old_gene_umi, bc_old_gene, bc_old_gene_umi)
                elif gene.endswith('--C'):
                    process_gene_data(barcode, geneid, umi, gene, new_gene, new_gene_umi, bc_new_gene, bc_new_gene_umi)
            except KeyError:
                continue

# Create expression matrices for new and old genes
def create_expression_matrix(gene_counts):
    """
    Create an expression matrix from the gene counts dictionary.
    """
    all_genes = list(set(gene.split('#')[0] for genes in gene_counts.values() for gene in genes))
    expression_matrix = pd.DataFrame(index=all_genes, columns=gene_counts.keys())

    for barcode, genes in gene_counts.items():
        for gene in genes:
            gene_name = gene.split('#')[0]
            expression_matrix.at[gene_name, barcode] = expression_matrix.get(gene_name, barcode, 0) + 1

    expression_matrix.fillna(0, inplace=True)
    return expression_matrix

new_gene_matrix = create_expression_matrix(bc_new_gene_umi)
old_gene_matrix = create_expression_matrix(bc_old_gene_umi)

new_gene_matrix.to_csv(f'{sample_folder}/new_gene_umi.csv')
old_gene_matrix.to_csv(f'{sample_folder}/old_gene_umi.csv')

combined_matrix = new_gene_matrix.add(old_gene_matrix, fill_value=0)
combined_matrix.to_csv(f'{sample_folder}/all_gene_umi.csv')

# Calculate UMI rates for new and old genes
def calculate_umi_rate(gene_umi_dict, cell_barcode):
    """
    Calculate the UMI rate for each barcode.
    """
    umi_rate = []
    for bc in cell_barcode:
        new_umi = len(gene_umi_dict.get(bc, {}))
        old_umi = len(bc_old_gene_umi.get(bc, {}))
        total_umi = new_umi + old_umi
        umi_rate.append(new_umi / total_umi if total_umi > 0 else 0)
    return umi_rate

# Calculate UMI rates for new and old genes
new_gene_umi_rate = calculate_umi_rate(bc_new_gene_umi, cell_barcode)
old_gene_umi_rate = calculate_umi_rate(bc_old_gene_umi, cell_barcode)

umi_rate_df = pd.DataFrame({
    'New Gene UMI Rate': new_gene_umi_rate,
    'Old Gene UMI Rate': old_gene_umi_rate
}, index=cell_barcode)

umi_rate_df.to_csv(f'{sample_folder}/gene_umi_rates_with_index.csv')

plt.figure(figsize=(4, 10))
box = plt.boxplot([new_gene_umi_rate, old_gene_umi_rate], patch_artist=True, medianprops={'color': 'black', 'linewidth': 2})
box['boxes'][0].set(facecolor='red', alpha=0.5)
box['boxes'][1].set(facecolor='blue', alpha=0.5)
plt.xticks([1, 2], ['New Gene', 'Old Gene'], fontsize=14)
plt.yticks(fontsize=14)
plt.ylim(0, 1)
plt.ylabel('Rate of Gene UMI Total Count per Cell', fontsize=14)

plt.savefig(os.path.join(sample_folder, "bc_new_old_gene_umi_rate_box_4.pdf"), bbox_inches='tight')
plt.savefig(os.path.join(sample_folder, "bc_new_old_gene_umi_rate_box_4.png"), bbox_inches='tight')
