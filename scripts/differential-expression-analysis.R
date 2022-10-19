# required libraries
library(VennDiagram)
library(edgeR)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(Rcpp)
library(EnhancedVolcano)


# required data
samples.metadata.clean <- readRDS(here::here("data","metadata","samples.metadata.RDS"))
filtered.counts <- readRDS(here::here("data","intermediate","filtered.counts.RDS"))

# Differential expression analysis
# Apply sample grouping based on sample_type from which the sample was derived

design <- model.matrix(~0+samples.metadata.clean$sample_type)
colnames(design) <- levels(as.factor(samples.metadata.clean$sample_type))

# colnames(design) <- c("tumor","normal")

# Estimate dispersions for tags
filtered.counts.dge <- estimateDisp(filtered.counts, design, robust = TRUE)

# Fit a generalized likelihood model to the DGELIST using sample grouping
fit <- glmFit(filtered.counts.dge,design)

#################################################################
# code in this section adapted from https://github.com/iscb-dc-rsg/2016-summer-workshop
# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(as.factor(samples.metadata.clean$sample_type)), 2))

comparisons <- list()
for (i in 1:nrow(condition_pairs)) {
  comparisons[[i]] <- as.character(condition_pairs[i,])
}

# remove MG to SG comparison
# comparisons[[2]] <- NULL

# vector to store deferentially expressed genes
sig_genes <- c()

# iterate over the contrasts, and perform a differential expression test for each pair
for (conds in comparisons) {
  # generate string contrast formula
  contrast_formula <- paste(conds, collapse=' - ') # should be (tumor - normal)
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=design)
  contrast_lrt <- glmLRT(fit, contrast=contrast_mat)
  topGenes <- topTags(contrast_lrt, n=Inf, p.value=0.05, adjust.method = "BH")
  
  # Grab highly ranked genes
  sig_genes <- union(sig_genes, rownames(topGenes$table))
}

# Filter out genes which were not differentially expressed for any contrast
de.genes <- filtered.counts.dge[rownames(filtered.counts.dge) %in% sig_genes,]
dim(de.genes$counts)
# 9697    5

################################################################

# Obtain the counts of genes expressed for each contrast individually
# This aims to obtain the number of genes differentially expressed between 
# the  i.e. tumor and normal

# tumor compared to normal
tumor_vs_normal <- glmLRT(fit, contrast=c(-1,1))


# Genes with most significant differences (using topTags)
topGenes <- topTags(tumor_vs_normal, adjust.method = "BH", p.value = 0.05, n=Inf)
dim(topGenes)
# 9697    5 

dim(topTags(glmLRT(fit, contrast=c(-1,1)), adjust.method = "BH", p.value = 0.05, n=Inf))
#######################################################################################
# DE genes at 5% FDR (using decideTestsDGE function)
#
# tumor compared to normal
tumor_vs_normal_de.genes <- decideTestsDGE(tumor_vs_normal, adjust.method = "BH", p.value = 0.05)

# get summary
summary(tumor_vs_normal_de.genes)
# -1*normal 1*tumor
# Down     3540
# NotSig   9295
# Up       6157

# create a dataframe with data on PV and SG differential gene expression
DE_data <- topGenes$table


DE_data <- tibble::rownames_to_column(DE_data, var = "gene_id")

write.csv(DE_data, here::here("results","tables","differentially_expressed_genes_all.csv"), row.names = FALSE)

# ensure to get only significant upregulated genes (pvalue < 0.05 and logFC < 2 and -2)
# get the significant up- and down-regulated genes
DE_data = subset(DE_data, abs(logFC) > 2)

# write out to excel upregulated and downregulated genes and the commmon genes between contrasts
DE_data <- DE_data[order(DE_data$logFC, decreasing = TRUE),]

write.csv(DE_data, here::here("results","tables","differentially_expressed_genes.csv"), row.names = FALSE)

########################################################################
# Plots
########################################################################
# Plotting to visually inspect differential gene expression results.

# Differential expression analysis - plots
#
# Volcano plots
tumor_vs_normal_DEGs = mutate(DE_data, sig=ifelse(DE_data$FDR <0.05 & logFC>2, "Up", ifelse(DE_data$logFC < -2, "Down", "Not Sig")))

p <- ggplot(tumor_vs_normal_DEGs, 
       aes(logFC, -log10(PValue))) +
  geom_point(aes(col=sig),size = 1) + 
  theme_bw(base_size = 9) + 
  #coord_cartesian(ylim=c(0,300))+ 
  coord_cartesian(xlim=c(-7,7)) +
  scale_color_manual(values=c("blue","black", "red")) +
  ggtitle("Tumor vs Normal differentially expressed genes.") 

p + geom_text_repel(data=filter(tumor_vs_normal_DEGs, logFC>4 | logFC< -5),
                  #family = "Times New Roman",
                  aes(label=gene_id),
                  #size = 2,
                  arrow = arrow(length = unit(0.02, 'npc')),
                  force = 7,box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.3, "lines"))

ggsave(here::here("results","figures","volcanoplot.png"), width = 7, height = 6)

#--------------------------------------------------------------
# Plot volcano plot using enhancedvolacano

top_upregulated_genes <- c("GRIA2","MYH2","PCDH10","CPB1","SCARNA22","CEACAM5","CYP2B7P")
top_downregulated_genes <- c("LEP", "GLYAT","PLIN4","TRARG1","HEPACAM","APOB","LOC107984268","LOC101926960")

selected_labels <- c(top_upregulated_genes, top_downregulated_genes)

keyvals.colour <- ifelse(DE_data$logFC < -2, 'skyblue', ifelse(DE_data$logFC > 2, 'red','black'))


EnhancedVolcano(DE_data,
                lab = rownames(DE_data),
                x = 'logFC',
                y = 'PValue',
                selectLab = selected_labels,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-32,
                FCcutoff = 2.0,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 1,
                colCustom = keyvals.colour,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
ggsave(here::here("results","figures","volcanoplot_enhanced.png"), width = 7, height = 6)

