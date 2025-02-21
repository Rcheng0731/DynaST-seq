# DynaST-seq
## Obtaining the Expression Matrix for DynaST-seq Data
1. Use the DynamicST tool (https://github.com/DynamicBiosystems/DynamicST) to process spatial transcriptomics data, converting raw image data and sequencing data into valuable gene expression information with a spatial context.
* The command for DynamicST count is:
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

