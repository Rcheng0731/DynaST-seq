# DynaST-seq
DynaST-seq: a Spatially and Temporally Resolved Transcriptomic Sequencing Approach for Profiling Gene Expression Dynamics in Tissues

## Abstract
Understanding the spatio-temporal gene expression is crucial for elucidating the molecular mechanisms that regulate cellular functions. Despite recent advances in spatial transcriptomics, most approaches lack temporal resolution, limiting insights into dynamic gene expression. Here, we present DynaST-seq, an advanced spatio-temporal transcriptomics technology that integrates metabolic RNA labeling with spatial transcriptomics, enabling high-resolution and transcriptome-wide profiling of RNA dynamics across space and time.

## Obtaining the Expression Matrix for DynaST-seq Data
1. Use the DynamicST tool (https://github.com/DynamicBiosystems/DynamicST) to process spatial transcriptomics data, converting raw image data and sequencing data into valuable gene expression information with a spatial context.
```
 DynamicST count \
    --sample sampleName \
    --id sampleName \
    --inputdir rawdata \
    --gtf xxx.gtf \
    --transcriptome xxx \
    --image HE.tif \
    --alignment alignment.json \
    --outputdir result
```
2. Consensus_sequence_construction.py：
   Used to group reads with identical UMIs, barcodes, and gene labels in the BAM file, and for each genomic position, select the most frequent base within each UMI group as the consensus base to construct the consensus sequence.
```
 python Consensus_sequence_construction.py \
    --bam_path xxx.bam \
    --gtf_file xxx.gtf \
    --temp_out_dir result
```
3. Expression_matrix_extraction.py：
   Extracts gene, barcode, and UMI information from the BAM file; calculates the UMI expression matrices and new-to-old ratios.
```
 python Expression_matrix_extraction.py \
    --sample_name sampleName \
    --sample_path result \
    --bam_file xxx.bam
```
4. Conversion_Proportions.py：
   Calculates and visualizations the substitution rates, and generates summary reports.
```
 python Conversion_Proportions.py \
    --sample_name sampleName \
    --gtf_file xxx.gtf \
    --bam_file xxx.bam \
    --out_dir result
```
5. Metabolic_Labeling_Correction.py：
   A binomial mixture model was used to correct the expression of new RNAs in spots of DynaST-seq data and estimate the substitution rate caused by metabolic RNA labeling.
```
 python Metabolic_Labeling_Correction.py \
    --mutations_file mutations_file \
    --reads_file reads_file \
    --group_info_file group_info_file \
    --convs_file convs.pkl.gz \
    --out_dir result
```

## Data
Raw data files are available at NCBI Gene Expression Omnibus (GEO)
