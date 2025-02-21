import os
import pysam
import ngs_tools as ngs
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

# Paths and configuration
sample_name = f"高低浓度对比/{sys.argv[1]}"
group_info = sys.argv[2]
bam_path = os.path.join('/home/disk/chengrui/workspace/时空转录组/data/可用数据汇总', sample_name, 'Aligned.sortedByCoord.UniqueGene.consensus.sorted.bam.tagST.bam')
out_dir = os.path.join('/home/disk/chengrui/workspace/时空转录组/data/可用数据汇总', sample_name)
temp_dir = os.path.join(out_dir, 'temp')
gtf_path = '/home/disk/chengrui/toolfolder/时空参考/mm10_gencodeM13/gencode.vM13.annotation.gtf'
n_threads = 34

barcodes = pd.read_csv(os.path.join(out_dir, 'outs/filtered_feature_bc_matrix/barcodes.tsv.gz'), header=None)[0].to_list()

os.makedirs(temp_dir, exist_ok=True)
gene_infos, transcript_infos = ngs.gtf.genes_and_transcripts_from_gtf(gtf_path, use_version=False)

split_bam_res = [value[0] for value in ngs.bam.split_bam(bam_path, temp_dir + "/temp", n=64, n_threads=n_threads, show_progress=True).values()]
for bam_file in split_bam_res:
    pysam.index(bam_file, f'{bam_file}.bai', '-@', str(n_threads))

def complement(base):
    """Returns the complementary base of a DNA nucleotide."""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return complement_map.get(base, base)

def process_bam_file(bam_file_path, barcodes):
    """Processes a BAM file and returns mutations and alignments data."""
    conversions, alignments, conversion_counts = [], [], defaultdict(int)
    with pysam.AlignmentFile(bam_file_path, 'rb', threads=n_threads) as bam_file:
        contigs = bam_file.references
        valid_contigs = {contig for contig in contigs if any(True for _ in bam_file.fetch(contig, until_eof=True))}
        
        for contig in valid_contigs:
            for record in bam_file.fetch(contig):
                barcode, gene_id, umi, strand = record.get_tag("CB").replace("_", ""), record.get_tag("GX"), record.get_tag("UB"), record.get_tag('ST')
                if barcode != "-" and gene_id != "-" and umi != "-" and barcode in barcodes:
                    base_counts, mutations = Counter(record.get_reference_sequence().upper()), {}
                    read_sequence = record.query_sequence.upper()
                    
                    for read_index, ref_index, ref_base in record.get_aligned_pairs(matches_only=True, with_seq=True):
                        read_base = read_sequence[read_index]
                        ref_base = ref_base.upper()

                        if 'N' in (ref_base, read_base):
                            continue

                        conversion = f'{ref_base}{read_base}' if strand == '+' else f'{complement(ref_base)}{complement(read_base)}'
                        conversion_counts[conversion] += 1

                        if ref_base != read_base:
                            mutations[ref_index] = {'read_id': record.query_name, 'contig': contig, 'position': ref_index, 'conversion': conversion}
                    
                    conversions.extend(mutations.values())
                    alignments.append({
                        'read_id': record.query_name, 'barcode': barcode, 'umi': umi, 'gene_id': gene_id,
                        'A': base_counts.get("A", 0), 'C': base_counts.get("C", 0),
                        'G': base_counts.get("G", 0), 'T': base_counts.get("T", 0),
                        'status': 'unassigned', 'assigned': bool(gene_id)
                    })
    
    return conversions, alignments, conversion_counts

def main(split_bam_res, barcodes, max_workers):
    """Main function to process all BAM files and gather results."""
    global all_conversions_list, all_alignments_list, conversion_counts
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_bam_file, bam_file_path, barcodes) for bam_file_path in split_bam_res]
        for future in futures:
            local_conversions, local_alignments, local_conversion_counts = future.result()
            all_conversions_list.extend(local_conversions)
            all_alignments_list.extend(local_alignments)
            for conversion, count in local_conversion_counts.items():
                conversion_counts[conversion] += count

    return all_conversions_list, all_alignments_list, conversion_counts

if __name__ == "__main__":
    max_workers = 40
    all_conversions_list, all_alignments_list, conversion_counts = main(split_bam_res, barcodes, max_workers)

    all_conversions_list.to_csv(os.path.join(out_dir, 'all_conversions_list.csv'), index=False)
    all_alignments_list.to_csv(os.path.join(out_dir, 'all_alignments_list.csv'), index=False)

    # Print conversion counts
    for conversion, count in conversion_counts.items():
        print(f"Conversion: {conversion}, Count: {count}")

    # Calculate total counts by start base
    start_base_totals = defaultdict(int)
    for conversion, count in conversion_counts.items():
        start_base_totals[conversion[0]] += count
    
    # Calculate conversion proportions
    conversion_proportions = {conversion: count / start_base_totals[conversion[0]] 
                            for conversion, count in conversion_counts.items() if conversion not in ['CC', 'AA', 'GG', 'TT']}
    
    labels = sorted(conversion_proportions.keys())
    proportions = [conversion_proportions[label] for label in labels]
    colors = plt.cm.rainbow(np.linspace(0, 1, len(labels)))

    plt.figure(figsize=(8, 4))
    plt.bar(labels, proportions, color=colors)
    plt.xlabel('Conversion Type')
    plt.ylabel('Proportion')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "conversion_proportions.pdf"), bbox_inches='tight')
    plt.show()

    pd.DataFrame(list(conversion_proportions.items()), columns=['Conversion_Type', 'Proportion']).to_csv(os.path.join(out_dir, 'conversion_proportions.csv'), index=False)
    pd.DataFrame(list(conversion_counts.items()), columns=['Conversion_Type', 'Count']).to_csv(os.path.join(out_dir, 'conversion_counts.csv'), index=False)
