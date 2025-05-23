library(tidyr)
library(dplyr)
library(Seurat)
library(cytosignal)
library(CellChat)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(future)
library(future.apply)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(org.Mm.eg.db)
library(ggrepel)
library(forcats)
library(aplot)

# Load CellChatDB mouse data
ccdb <- CellChatDB.mouse
def_g2u <- ccdb$geneInfo[, c('Symbol', 'EntrezGene.ID')] %>%
  na.omit() %>%
  rename(gene_name = Symbol, uniprot = EntrezGene.ID)

# Subset CellChatDB interaction data
def.index <- ccdb$interaction
lig.df <- def.index[, c('interaction_name_2', 'ligand', 'annotation')]
rec.df <- def.index[, c('interaction_name_2', 'receptor', 'annotation')]

# Process ligand and receptor information
process_ligand_receptor <- function(df, annotation_type) {
  df %>%
    filter(annotation %in% annotation_type) %>%
    separate_rows(ligand, sep = ", ") %>%
    distinct()
}

lig.df <- process_ligand_receptor(lig.df, c('Secreted Signaling', 'Non-protein Signaling'))
rec.df <- process_ligand_receptor(rec.df, c('Secreted Signaling', 'Non-protein Signaling'))

# Create ligands and receptors factor with appropriate levels
create_factor <- function(df, annotation_type, col_name) {
  df %>%
    filter(annotation %in% annotation_type) %>%
    mutate(
      name = factor(.[[col_name]]),
      ligand_or_receptor = factor(.[[col_name]], levels = unique(.[[col_name]]))
    ) %>%
    select(interaction_name_2, ligand_or_receptor)
}

lig <- create_factor(lig.df, c('Secreted Signaling', 'Non-protein Signaling'), 'ligand')
rec <- create_factor(rec.df, c('Secreted Signaling', 'Non-protein Signaling'), 'receptor')

# Combine ligand and receptor info into a list
combine_lig_rec <- function(lig, rec) {
  list(
    combined = c(lig, rec),
    ligands = lig,
    receptors = rec
  )
}

defdb.diff <- combine_lig_rec(lig, rec)

# Prepare for contact type
def.index_no_underscore <- def.index %>%
  filter(!grepl("_", ligand) & !grepl("_", receptor)) %>%
  filter(!grepl(" ", ligand) & !grepl(" ", receptor))

# Create ligand-receptor pair list
ligand_receptor_pairs <- mapply(function(ligand, receptor) c(ligand, receptor), 
                                def.index_no_underscore$ligand, 
                                def.index_no_underscore$receptor, 
                                SIMPLIFY = FALSE)

write.csv(def.index_no_underscore, "/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/4_MCC/ref_data/1_cellchatdb_1v1.csv")

region <- "Inner_stripe"
adata_sec <- readRDS(paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/data/4_MCC/", region, ".rds"))
metadata <- read.csv(paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/4_MCC/total_region_res/", region, "/coloc_results_meta.csv"), row.names = 1)
label <- "all_lr"

process_sample_data <- function(adata_sec_sub, sample, LRpair, knn = 8, normalize = TRUE) {
  newdata <- GetAssayData(adata_sec_sub, assay = "NewRNA", slot = "counts")
  olddata <- GetAssayData(adata_sec_sub, assay = "OldRNA", slot = "counts")
  count <- newdata + olddata
  location <- as.data.frame(adata_sec_sub@images[[sample]]$centroids@coords[, c('x', 'y')])
  
  countsum <- Matrix::colSums(count)
  expr <- count
  if (normalize) {
    expr <- Matrix::t(log(Matrix::t(expr)/countsum * median(countsum) + 1))
  }
  ligand <- expr[LRpair[1], ]
  receptor <- expr[LRpair[2], ]
  LRexp <- rbind(ligand, receptor)
  nnmatrix <- RANN::nn2(location, k = knn)$nn.idx
  
  neighexp <- apply(nnmatrix, 1, function(x) apply(LRexp[, x[2:knn]], 1, max))
  LRadd <- pmax(LRexp[1, ] * neighexp[2, ], LRexp[2, ] * neighexp[1, ])
  
  LRadd_max <- quantile(LRadd, probs = 0.95)
  LRadd[LRadd > LRadd_max] <- LRadd_max
  alpha <- (LRadd - min(LRadd)) / (max(LRadd) - min(LRadd)) * (1 - 0.1) + 0.1
  
  list(LRadd = LRadd, LRexp = LRexp, location = location, alpha = alpha)
}

all_tmp_list <- list()
all_tmpLRadd_list <- list()

pdf(paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/5_spagene/", region, "/", label, "/medium_lr.pdf"), width = 35, height = 8)
for (lr in ref$lr) {
  print(lr)
  plots_list <- list()
  parts <- strsplit(lr, "_")[[1]]
  l <- parts[1]
  r <- parts[2]
  LRpair <- c(l, r)

  for (sample in c('control_1', 'control_2', 'IRI_4h_1', 'IRI_4h_2', 'IRI_12h_1', 'IRI_12h_2', 'IRI_24h_1', 'IRI_24h_2')) {
    adata_sec_sub <- subset(adata_sec, orig.ident == sample)
    processed_data <- process_sample_data(adata_sec_sub, sample, LRpair)
    
    tmp <- data.frame(x = processed_data$location[, 1], 
                      y = processed_data$location[, 2], 
                      Exp = as.factor(expcol), 
                      sample = sample, lr = lr)
    
    all_tmp_list[[paste(sample, lr, sep = "_")]] <- tmp
    tmpLRadd <- data.frame(x = processed_data$location[, 1], 
                           y = processed_data$location[, 2], 
                           LR = processed_data$LRadd, 
                           sample = sample, lr = lr)
    
    all_tmpLRadd_list[[paste(sample, lr, sep = "_")]] <- tmpLRadd
  }
  
  print(plots_list[[1]] / plots_list[[5]] | plots_list[[2]] / plots_list[[6]] | plots_list[[3]] / plots_list[[7]] | plots_list[[4]] / plots_list[[8]])
}
dev.off()

all_tmp_df <- do.call(rbind, all_tmp_list)
all_tmpLRadd_list_df <- do.call(rbind, all_tmpLRadd_list)

write.csv(all_tmp_df, paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/5_spagene/", region, "/", label, "/medium_all_tmp_df.csv"), row.names = FALSE)
write.csv(all_tmpLRadd_list_df, paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/5_spagene/", region, "/", label, "/medium_all_tmpLRadd_list_df.csv"), row.names = FALSE)

all_tmpLRadd_list_df_res <- all_tmpLRadd_list_df %>%
  group_by(lr, sample) %>%
  summarize(avg_LR = mean(LR), .groups = 'drop')

all_tmpLRadd_list_df_res_wide_df <- all_tmpLRadd_list_df_res %>%
  pivot_wider(names_from = sample, values_from = avg_LR, values_fill = list(avg_LR = 0))

all_tmpLRadd_list_df_res_wide <- data.frame(
  control = rowMeans(all_tmpLRadd_list_df_res_wide_df[, c("control_1", "control_2")]),
  IRI_4h = rowMeans(all_tmpLRadd_list_df_res_wide_df[, c("IRI_4h_1", "IRI_4h_2")]),
  IRI_12h = rowMeans(all_tmpLRadd_list_df_res_wide_df[, c("IRI_12h_1", "IRI_12h_2")]),
  IRI_24h = rowMeans(all_tmpLRadd_list_df_res_wide_df[, c("IRI_24h_1", "IRI_24h_2")])
)

pdf(paste0("/home/disk/chengrui/workspace/时空转录组/DynamicST/最终需求的汇总脚本/9_细胞通讯/result/5_spagene/", region, "/", label, "/medium_lr_score_scale.pdf"), width = 5, height = 7)
print(pheatmap(t(scale(t(all_tmpLRadd_list_df_res_wide))),
               scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, 
               show_rownames = TRUE, show_colnames = TRUE, 
               display_numbers = TRUE))
dev.off()