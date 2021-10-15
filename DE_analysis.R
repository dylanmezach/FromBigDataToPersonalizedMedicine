#### Differential gene expression analysis ###
# This script will perform a differential expression analysis of 2 different conditions.
# It will also produce a manhattenplot, highlighting genes targeted by a transcription
# factor of wich the motif is encriched in the DE genese 

### Required packages ###
library(DESeq2)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(biomaRt)
library(qqman)
library(dplyr)
library(GLAD)

# Import data #
mRNA_data <- read.table('g:/Mijn Drive/FromBigDataToPersonalizedMedicine/Data/BMS_practicum_okt2021.expression.genelevel.v75.htseq.txt', header = T, row.names = 1)

# Use the information in the column names to generate a data.frame comprising the sample information.
colData <- data.frame(row.names = colnames(mRNA_data), 
                      sample= unlist(lapply(colnames(mRNA_data), 
                                            function(x){unlist(strsplit(x, split = "_"))[1]})), 
                      stimulation= unlist(lapply(colnames(mRNA_data), 
                                                 function(x){unlist(strsplit(x, split = "_"))[2]})))

# Add library size, which is the total ammount of gene reads per sample
colData$libSize <- colSums(mRNA_data)
#First 10 rowns of the colData object
head(colData)

####### PCA #######

dds <- DESeqDataSetFromMatrix(countData = mRNA_data,
                              colData = colData,
                              design = ~ sample)

vst <- assay(vst(dds, blind=FALSE))
# To calculate the components by sample we need to transpose our matrix of normalized gene expression 
pcData <- prcomp(t(vst))
pcVar <- summary(pcData)

varPC1 <- pcVar$importance[2,1]
varPC2 <- pcVar$importance[2,2]
pcPlotData <- data.frame(pcData$x[,1:4], colData[rownames(pcData$x),])

PCA_group <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=sample))+
  geom_jitter(alpha=0.6, aes(size = 1))+
  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
  scale_color_nejm()+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Variation between groups")+
  guides(col = guide_legend(ncol = 8))

PCA_stim <- ggplot(pcPlotData, aes(x=PC1 , y=PC2 , color=stimulation))+
  geom_jitter(alpha=0.6, aes(size=1))+
  xlab(paste0("PC1 explained variance = ", varPC1*100, "%"))+
  ylab(paste0("PC2 explained variance = ", varPC2*100, "%"))+
  scale_color_nejm()+
  theme_bw()+
  theme(legend.position = "bottom")+
  ggtitle("Variation between stimulations")+
  guides(col = guide_legend(ncol = 8))


tiff('G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/PCA.tiff', units = "in", height = 5, width = 15, res = 300)
grid.arrange(PCA_group, PCA_stim, nrow=1)
dev.off()


#select the sample names of 2 different stimulations
colData_stim2_4_names <- rownames(colData)[which(colData$stimulation == 'redBDOctrl' | colData$stimulation == 'redBDOIL15IL21IFNb')]
#select the read counts of only the 2 selected stimulations
stim2_4_countTable <- mRNA_data[, colData_stim2_4_names]
# make a colData with only the 2 selected stimulations
stim2_4_colData <- colData[colData_stim2_4_names,]

# Make input data into right format for DESeq2 
dds_stim_2_4 <- DESeqDataSetFromMatrix(countData = stim2_4_countTable, colData = stim2_4_colData, design = ~ stimulation)
# perform normalization
vst_stim_2_4 <- assay(vst(dds_stim_2_4, blind = F))
# Perform DE analyis
dds_DE_stim_2_4 <- DESeq(dds_stim_2_4)
# Extract results from DE analysis
DE_results <- as.data.frame(results(dds_DE_stim_2_4))

####### Create control genese for TFBSA ####
sub_res <- rownames(DE_results[order(DE_results$padj, decreasing = T),])[1:69]
write.csv(sub_res, "G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Data/control_Genes.csv")
# save DE results
write.csv(DE_results[,c(2,6)], "G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/DE_results_all_genes.csv")


# Make stimulation column into factor data
stim2_4_colData$stimulation <- factor(stim2_4_colData$stimulation)
# select DE genes with a adjusted-P-value of less then 0.05
DE_stim_2_4 <- DE_results[which(DE_results$padj <= 0.05),]
# show number of DE genes
nrow(DE_stim_2_4)

# Sort genes based on adjusted-p-val
sorted_DE_stim_2_4 <- DE_stim_2_4[order(DE_stim_2_4$padj),]
# make dataframe of only genecodes, log2FC and adjp
sub_sorted_DE_stim_2_4 <- as.data.frame(sorted_DE_stim_2_4[,c("log2FoldChange", "padj")])
# save to csv file in Gdrive
write.csv(sub_sorted_DE_stim_2_4, 'G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/sorted_DE_stim2_4.csv')



### creating manhatten plot ####
#use Ref genome annotation info from biomart
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
# select the attributes you want
genemap <- as.data.frame(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(DE_results),
                 mart = ensembl ))

# remove duplicate rows
genemap <- genemap[!duplicated(genemap$ensembl_gene_id),]

# make ensemble codes to rownames
rownames(genemap) <- genemap$ensembl_gene_id
genemap$ensembl_gene_id <- NULL

# make list of significant DE genes
sig_names <- rownames(sub_sorted_DE_stim_2_4)

# merge gene map with DE results to comine padj and log2FC
genemap <- merge(genemap, DE_results, by= 0)
# convert character chromosome names to numeric chromosome indicators (x = 23, y = 24)
genemap$chromosome_name <- ChrNumeric(genemap$chromosome_name)
#removal of NA's
genemap <- genemap[complete.cases(genemap),]
#subset top 7 DE genes

# make vulacano plot to visualize DE genes
tiff("G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/VulcanoPlot.tiff", units = "in", height = 5, width = 5, res = 300)
ggplot(genemap, aes(x=log2FoldChange, y= -log10(padj)))+
  geom_point(aes(color= padj <= 0.05))+
  geom_hline(yintercept = 0, lwd=1, alpha= 0.6)+
  geom_vline(xintercept = 0, lwd=1, alpha= 0.6)+
  scale_color_d3()+
  ggtitle("DE between redBDO and redBDO+IL-15+IL-21+IFNb")+
  theme_bw()
dev.off()

#read in results for transcription factor binding site
TFBSA_results <- read.table("G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Data/TF binding site/fimo.tsv", sep =  "\t", header = T)
# read in results of differential expression analysis
DE_results <- read.table("G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/DE_results_all_genes.csv", row.names = 1, sep = ',', header = T)
#create a genemap with needed attributes
genemap <- as.data.frame(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                               filters = "ensembl_gene_id",
                               values = rownames(DE_results),
                               mart = ensembl ))
# separate the chr no, start pos and end pos of TFBSA
motif_pos <- as.data.frame(t(as.data.frame(apply(TFBSA_results, 1, function(motif){
  tmp <- strsplit(motif[3], split = ":")
  tmp2 <- strsplit(tmp$sequence_name[2], split = "-")
  return(c(tmp$sequence_name[1], tmp2[[1]]))
}))))

comb <- cbind(TFBSA_results,motif_pos)
colnames(comb)[11:13] <- c("chromosome_name", "start_pos", "end_pos")
# convert chr names to numeric values
comb$chromosome_name <- ChrNumeric(comb$chromosome_name)

#add 350 bp to the start site to make it the start site of the gene instead of promotor region
comb$start_pos <- as.integer(comb$start_pos)+350
# merge TFBSA results with genemap info to retrieve genes of interest(GOI)
test <- merge(genemap, comb, by.x = 'start_position', by.y = 'start_pos')
#select only unique genes of interest
GOI <- unique(test$hgnc_symbol)

target <- merge(DE_results, genemap, by.x = 0, by.y = "ensembl_gene_id")
target$chromosome_name <- ChrNumeric(target$chromosome_name)
target <- target[complete.cases(target),]
# create manhatten plot with highlighting genes that were enriched in the TFBSA and write to Gdrive
tiff("G:/Mijn Drive/FromBigDataToPersonalizedMedicine/Results/ManhattanPlot.tiff", units = 'in', height = 5, width = 11, res = 300 )
manhattan(target,
          chr = "chromosome_name",
          bp = 'start_position', 
          p= "padj", 
          snp = 'hgnc_symbol', 
          suggestiveline = F, 
          genomewideline = F,
          highlight = GOI,
          annotatePval = T,
          main = "DE genes arranged on chromosomal position")
dev.off()


