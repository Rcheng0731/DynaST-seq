# DynaST-seq
DynaST-seq: a Spatially and Temporally Resolved Transcriptomic Sequencing Approach for Profiling Gene Expression Dynamics in Tissues
## Obtaining the Expression Matrix for DynaST-seq Data
1. Use the DynamicST tool (https://github.com/DynamicBiosystems/DynamicST) to process spatial transcriptomics data, converting raw image data and sequencing data into valuable gene expression information with a spatial context.
```
 DynamicST count \
    --sample sampleName \
    --id sampleName \
    --inputdir rawdata \
    --gtf Homo_sapiens.GRCh38.99.gtf \
    --transcriptome Homo_sapiens_GRCh38 \
    --image HE.tif \
    --alignment alignment.json \
    --outputdir result
```
2. Consensus_sequence_construction.py：used to retain reads with identical UMIs, barcodes, and gene labels in the BAM file, and for each genomic position, select the most frequent base within each UMI group as the consensus base to construct the consensus sequence.
3. Expression_matrix_extraction.py：extracts gene, barcode, and UMI information from the BAM file, calculates the UMI expression matrices and UMI rates for new and old genes.
4. Conversion_Proportions.py：calculates mutation and conversion counts, and generates summary reports and visualizations of conversion proportions for sequencing data analysis.
5. Metabolic_Labeling_Correction.py：binomial mixture model was used to correct the distribution of T-to-C substitutions in spots of DynaST-seq data and estimate the substitution rate caused by metabolic RNA labeling.
