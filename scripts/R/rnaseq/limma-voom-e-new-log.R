# Perform differential expression from limma-voom (with eBayes) 

rm(list = ls())

load("readcounts_dfs.RData")

library(edgeR)

# use edgeR to remove lowly expressed genes and normalize reads for
# sequencing depth ; see code chunks above
info_e.edger <- factor(c(rep("4C", 3) , rep("4C11+", 3)))
info_e.edger <- relevel(info_e.edger , ref = "4C")
edgeR.DGElist <- DGEList(counts = df_e , group = info_e.edger)
cpm <- cpm(edgeR.DGElist) # log took off because of Diogo's script

keep <- rowSums(cpm > 1) >= 0.35 * length(colnames(counts)) # this was changed according to the script of Diogo Pessoa on August 20th 2018

edgeR.DGElist <- edgeR.DGElist[keep ,]
edgeR.DGElist <- calcNormFactors(edgeR.DGElist, method = "TMM")

# limma also needs a design matrix , just like edgeR
design <- model.matrix(~info_e.edger)

# transform the count data to log2 -counts -per - million and estimate
# the mean - variance relationship , which is used to compute weights
# for each count -- this is supposed to make the read counts
# amenable to be used with linear models
rownames(design) <- colnames(edgeR.DGElist)
voomTransformed <- voom(edgeR.DGElist , design, plot = FALSE )

# fit a linear model for each gene
voomed.fitted <- lmFit(voomTransformed, design = design)

# compute moderated t- statistics , moderated F- statistics ,
# and log - odds of differential expression
voomed.fitted <- eBayes(voomed.fitted )

# extract gene list with logFC and statistical measures
colnames(design) # check how the coefficient is named
DGE.results_limma <- topTable(voomed.fitted, coef = "info_e.edger4C11+",
                                   number = Inf , adjust.method = "BH", sort.by = "logFC")

# remove the dot containing the version of ENSEMBL GENE (with the dot the gene converter does not recognize the ENSEMBL gene id)
row.names(DGE.results_limma) <- gsub(pattern = "\\..{1,}", "", row.names(DGE.results_limma)) 

# get the gene symbol annotation of the ENSEMBLGENE genes
library(Mus.musculus)

# separate only the ENSEMBLGENE
DGEgenes <- rownames(DGE.results_limma)

# get gene symbol using mapIds, because with this it get only uniq IDs
anno <- select(org.Mm.eg.db, keys = DGEgenes, keytype = "ENSEMBL", columns = c("SYMBOL", "GENENAME"))

# transform anno in data frame
anno.df <- as.data.frame(anno)

# bind anno.df and DGE.results_limma by ENSEMBLGENE using "merge" = more accurate than cbind
DGE.results_limma_anno <- merge(anno.df, DGE.results_limma, by.x = "ENSEMBL",by.y = "row.names")
# order according to logFC (for the preranked list)
DGE.results_limma_anno <- DGE.results_limma_anno[order(DGE.results_limma_anno$logFC), ]

# save one version with NAs
write.table(DGE.results_limma_anno, "DGE.results_limma_anno_e.txt", quote = F,
              row.names = F, sep = "\t")

# remove NAs
# I installed the package rgr to get this function
DGE.results_notNA <- rgr::remove.na(DGE.results_limma_anno)
DGE.results_notNA <- DGE.results_notNA$x
DGE.results.sorted_notNA <- DGE.results_notNA[order(DGE.results_notNA$logFC), ]
# save version without NAs
write.table(DGE.results_notNA, "DGE.results_limma_anno_notNA_e.txt", 
            quote = F,
            row.names = F, sep = "\t")

# save subset only with gene symbols and logFC (for GSEA preranked)
preranked <- subset(DGE.results_notNA, select= c("SYMBOL", "logFC"))
write.table(preranked, "preranked_geneexpr_e.txt", quote = F,
            row.names = F, sep = "\t")


save.image("limma-voom-e.RData")

# venn diagram for the pairwise comparisons
# DE_list <- list(a = rownames(subset(DGE.results_limma , adj.P.Val <= 0.05)),
#                 b = rownames(subset(DGE.results_limma , adj.P.Val <= 0.05)),
#                 c = rownames(subset(DGE.results_limma , adj.P.Val <= 0.05)))
#                  
# gplots::venn(DE_list)
# 
# DE_gns <- UpSetR :: fromList(DE_list)
# UpSetR::upset (DE_gns , order.by = "freq")
