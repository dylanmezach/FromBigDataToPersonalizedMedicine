#====From Big Data to Personalized Medicine====

#----PATHWAY ANALYSIS----
#----1. PREPARATION----
# These packages are from CRAN
install.packages('data.table')
install.packages('ggplot2')
install.packages('dplyr')
install.packages("ggsci")
install.packages("ggrepel")

# These packages are from Bioconductor
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
BiocManager::install("topGO")
BiocManager::install("DOSE", force=TRUE)
BiocManager::install("DESeq2")

# Library
library(data.table)
library(ggplot2)
library(dplyr)
library(ggsci)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(topGO)
library(DOSE)
library(DESeq2)

set.seed(123)

#Import dataset
library(readxl)
Book1 <- read.table("/My Downloads/DE_results_all_genes.csv", sep = ",", header = T)
View(Book1)
summary(Book1)

#----2. PREPARE THE DATA----
# convert HGNC gene symbols to more stable ENTREZ IDs
entrez <- bitr(Book1$X, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
head(entrez)
dim(entrez)

# Convert ENTREZ to text, so that R would process it correctly
entrez$ENTREZID <- as.character(entrez$ENTREZID)

# Add column with ENTREZ IDs to the table, remove rows where ENTREZ was missing
Book1 <- merge(Book1, entrez, by.x = 'X', by.y = 'ENSEMBL')
View(Book1)
dim(Book1)

# Prepare the data for overrepresentation tests
# Filter in significant (adjusted P<0.05) results
# Filter in only genes with larger than than 2-fold expression change 
Book1_sig <- Book1[which(Book1$padj < 0.05 & abs(Book1$log2FoldChange) > 1), ]
View(Book1_sig)

# Order gene table by effect size (log2 fold-change), from largest to smallest
Book1 <- Book1[order(Book1$log2FoldChange, decreasing = T), ]

# Look at the range of fold changes
plot(Book1$log2FoldChange, xlab = 'Gene', ylab = 'log2(FC)')

# Prepare data for GSEA (vector of log2(FC)'s, named by ENTREZ IDs)
Book1_gsea <- Book1$log2FoldChange
names(Book1_gsea) <- Book1$ENTREZID

Book1_gsea <- sort(Book1_gsea, decreasing = T)
head(Book1_gsea)
# This vector of gene names corresponds to all genes we tested in the analysis and will be used later as "gene universe"


#----3. OVER REPRESENTATION ANALYSIS (ORA)----

#----3.1 RUN KEGG OVER-REPRESENTATION ANALYSIS----
#In this part we run enrichment test (hypergeometric test) for all differentially expressed genes (FDR<0.05, FC>2).
#We use all genes tested in the RNA-profiling study as a background set or "gene universe". This list was already constructed in the previous step.
#We will query for all enrichment results and write out 25 most differentially expressed, not accounting the significance.
#The default multiple testing is done by Benjamini-Hochberg method, flag pvalueCutoff can be used for filtering only significant results (<0.05). Additionally, qvalueCutoff flag is used to filter results based on another popular method of FDR estimation (Storey q-value).
#The defaults of the command also apply restrictions on the sizes of gene sets tested in the analysis, those can be seen with command ?enrichKEGG.

KEGG_all <- enrichKEGG(gene = Book1_sig$ENTREZID,
                       organism = 'hsa',
                       universe = Book1$ENTREZID,
                       pvalueCutoff = 1, 
                       qvalueCutoff = 1)
KEGG_all <- setReadable(KEGG_all, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
View(KEGG_all@result)

View(KEGG_all@result)
write.table(KEGG_all@result[1:17,], 'KEGG_OR_Table.csv', sep = ',', row.names = F)


# show the data structure
head(str(KEGG_all),25)

# result table
head(KEGG_all@result)  
# what ontology was tested
KEGG_all@ontology  
# what P-value cutoff was used
KEGG_all@pvalueCutoff     

#Now we look at the results by printing out 25 first rows of the result table. 
View(head(KEGG_all@result, 17))
View(head(KEGG_all@result, 5))
View(head(KEGG_all@result), )


#----3.1.1 VISUALIZE OVERALL KEGG ENRICHMENT ANALYSIS RESULTS----
input_barplot <- KEGG_all@result
input_barplot$Description <- factor(input_barplot$Description, levels = rev(as.character(input_barplot$Description)))
# here we apply the default significance thresholds (Benjamini-Hochberg P<0.05 and Storey Q-value<0.2)
input_barplot <- input_barplot[input_barplot$p.adjust < 0.05 & input_barplot$qvalue < 0.2, ]

BP17 <- ggplot(input_barplot, aes(x = Description, y = -log10(pvalue), fill = Count)) + geom_bar(stat = 'identity') + 
  theme_classic() + 
  coord_flip() + scale_fill_continuous(low = 'lightblue', high = 'salmon') +
  ggtitle("KEGG over-representation analysis")
BP17

tiff("./KEGG_barplot_17_cat.tiff", units = "in", height = 10, width = 10, res = 300)
BP17
dev.off()

#Next we visualize five most enriched pathways on the network graph. 
#Size of the pathway node shows statistical significance and links indicate gene membership in the pathway. 
#This gives the static graph as a output. 
#If it is necessary to investigate further, you can use flag fixed = FALSE to see interactive version.

png('~/R/KEGG_cnetplot.png', width = 20, height = 20, res = 400, units = 'in')
Networkgraph <- cnetplot(KEGG_all, categorySize = "pvalue", showCategory = 5, foldChange = Book1_gsea)
dev.off()
Networkgraph

# Interactive version
png('~/R/KEGG_cnetplot2.png', width = 20, height = 20, res = 400, units = 'in')
cnetplot(KEGG_all, categorySize = "pvalue", showCategory = 5, foldChange = Book1_gsea, fixed = F)
dev.off()


#----3.1.2 VISUALIZE THE MOST SIGNIFICANTLY ENRICHED KEGG PATHWAY
#For KEGG pathways it is possible to visualize the location of differentially expressed genes on the pathway, as well as their magnitude of differential expression. 
#This kind of visualization might help to identify most relevant genes for given condition and generate new research hypotheses.
# Following command saves the corresponding KEGG graph to your work directory

#For most significant pathway - Influenza A
pathview(gene.data = Book1_gsea,
         pathway.id = KEGG_all@result[1, ]$ID,
         species = "hsa",
         limit = list(gene = max(abs(Book1_gsea)), cpd = 1), 
         out.suffix = "most_sig_pathway")

#For NOD-like receptor signaling pathway
pathview(gene.data = Book1_gsea,
         pathway.id = KEGG_all@result[4, ]$ID,
         species = "hsa",
         limit = list(gene = max(abs(Book1_gsea)), cpd = 1), 
         out.suffix = "most_sig_pathway")


#----3.2 RUN GENE ONTOLOGY OVER-REPRESENTATION ANALYSIS----
#Next we run the overrepresentation test for Gene Ontologies.
#For the sake of brevity, we test only the ontologies from Biological Process category. 
#By setting the ont flag, is also possible to test Molecular Function (MF), Cellular Component (CC) or all three categories combined (ALL).
#For visualization purposes we apply this time significance thresholds (Benjamini-Hochberg FDR<0.05 and Storey FDR<0.05).
#We specify that we test the enrichment for the gene sets with sizes between 10-10000 genes (flags minGSSize and maxGSSize)

#GENE-IDs
GO_all_ID <- enrichGO(gene = Book1_sig$ENTREZID,
                   universe = Book1$ENTREZID,
                   OrgDb = 'org.Hs.eg.db',
                   ont = 'BP',
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05,
                   minGSSize = 10, 
                   maxGSSize = 10000)

GO_all_ID <- setReadable(GO_all_ID, OrgDb = org.Hs.eg.db)

View(GO_all_ID@result)
View(head(GO_all_ID@result, 17))

GO_all_ID_table <- data.frame(GO_all_ID[1:17,])
write.table(GO_all_ID_table, file = "~/R/GO_all_table.txt", sep = ",", quote = FALSE, row.names = F)

#Dotplot
dotplot(GO_all_ID) + ggtitle("Dotplot")
no_sig <- nrow(GO_all_ID@result[which(GO_all_ID@result$p.adjust < 0.0001 & GO_all_ID@result$qvalue < 0.2),])
#GENE-SYMBOL
GO_all_symbol <- enrichGO(gene = Book1_sig$ENTREZID,
                   universe = Book1$ENTREZID,
                   OrgDb = 'org.Hs.eg.db',
                   ont = 'BP',
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.05,
                   minGSSize = 10, 
                   maxGSSize = 10000, 
                   readable = TRUE)

View(head(GO_all_symbol@result))

write.table(GO_all_symbol@result[1:17,], 'GO_OR_Table.csv', sep = ',', row.names = F)

#Dotplot
dotplot(GO_all_symbol, showCategory=no_sig) + ggtitle("Dotplot")

#Heatmap
heatplot(GO_all_symbol, foldChange=Book1_gsea) + ggtitle("Heatmap")

#Barplot
input_barplot_GO <- GO_all_ID@result[1:17,]
input_barplot_GO$Description <- factor(input_barplot_GO$Description, levels = rev(as.character(input_barplot_GO$Description)))
# here we apply the default significance thresholds (Benjamini-Hochberg P<0.05 and Storey Q-value<0.2)
input_barplot_GO <- input_barplot_GO[order(input_barplot_GO$p.adjust, decreasing = T),]
input_barplot_GO <- input_barplot_GO[input_barplot_GO$p.adjust < 0.05 & input_barplot_GO$qvalue < 0.2, ]
GO17 <- ggplot(input_barplot_GO, aes(x = Description, y = -log10(pvalue), fill = Count)) + geom_bar(stat = 'identity') + 
  theme_classic() + 
  coord_flip() + scale_fill_continuous(low = 'lightblue', high = 'salmon') +
  ggtitle("GO over-representation analysis")


tiff("./GO_barplot_17_cat.tiff", units = "in", height = 8, width = 8, res = 300)
GO17
dev.off()

#----3.2.1 VISUALIZE THE GO ENRICHMENT ANALYSIS RESULTS----
#To ease the interpretation and see the hierarchical relationships between most enriched GO terms, we visualize those on graph structure. 
#Square boxes are significantly enriched terms and color depicts the significance (P-value). 
#We will visualize 15 most significant terms.

#GENE-ID
png('~/R/GO_plot_ID.png', height = 10, width = 10, units = 'in', res = 400)
plotGOgraph(GO_all_ID, firstSigNodes = 5, sigForAll = TRUE)
dev.off()

#GENE-SYMBOL
png('~/R/GO_plot_symbol.png', height = 10, width = 10, units = 'in', res = 400)
plotGOgraph(GO_all_symbol, firstSigNodes = 5, sigForAll = TRUE)
dev.off()


#----3.3 RUN DISEASE ONTOLOGY OVER-REPRESENTATION ANALYSIS----
#Clusterprofiler allows to do also ORA for human diseases, making use of several databases of curated gene-disease relationships. 
#As our trait of interest is Coeliac disease, it would be interesting to see differentially expressed genes also show enrichment of some human diseases and whether those are could be related with Coeliac disease or any other autoimmune disease.

DO_all <- enrichDO(gene = Book1_sig$ENTREZID,
                   universe = Book1$ENTREZID,
                   pvalueCutoff = 1, 
                   qvalueCutoff = 1,
                   ont = "DO")
View(head(DO_all@result, 25))


#----4. GENE SET ENRICHMENT ANALYSIS (GSEA)----
#Next we will use GSEA (Subramanian et al. 2005) on the differential analysis results. 
#Unlike ORA, GSEA uses all the results from differential expression analysis, not only significant ones.

#----4.1 RUN GSEA for KEGG pathways
#For the sake of speed and brevity, we use recommended minimal number of 1000 permutations for this analysis.
#To get more stable and precise results, your could increase the number of permutations to e.g. 10000 in real scientific work. 
#Also, we test only KEGG pathways which have more than 50 known gene members.

#GSEA_KEGG
KEGG_GSEA <- gseKEGG(geneList = Book1_gsea,
                     organism = 'hsa',
                     nPerm = 10000,
                     minGSSize = 50,
                     pvalueCutoff = 1,
                     verbose = FALSE)

View(head(KEGG_GSEA@result, 1000))
View(KEGG_GSEA@result)

write.table(KEGG_GSEA@result[1:5,], 'KEGG_GSEA_Table.csv', sep = ',', row.names = F)


#GSEA_GO
GO_GSEA <- gseGO(geneList = Book1_gsea,
                 ont = 'BP', 
                 OrgDb = org.Hs.eg.db,
                 nPerm = 10000,
                 keyType = 'ENTREZID',
                 minGSSize = 50,
                 pvalueCutoff = 1,
                 verbose = FALSE)

View(GO_GSEA@result)

write.table(GO_GSEA@result[1:5,], 'GO_GSEA_Table.csv', sep = ',', row.names = F)


#----4.1.1 VISUALIZE GSEA PLOT FOR MOST SIGNIFICANT KEGG PATHWAY----

png('~/R/GSEA_plot.png', height = 6, width = 9, units = 'in', res = 400)

#Here the geneSetID we specify the first (most significant) GSEA result to visualize:
gseaplot(KEGG_GSEA, geneSetID = head(KEGG_GSEA@result, 1)$ID)
dev.off()

#GSEA result of NOD-like receptor signaling pathway
png('~/R/GSEA_plot_NOD.png', height = 6, width = 9, units = 'in', res = 400)
gseaplot(KEGG_GSEA, geneSetID = "hsa04621", title = "GSEA NOD-like receptor signaling pathway")
dev.off()

GSEA_NOD <- gseaplot(KEGG_GSEA, geneSetID = "hsa04621")
GSEA_NOD

#leading edge analysis
library(DOSE)

#Visualize first 5 KEGG pathways in cnetplot
KEGG_GSEA2 <- setReadable(KEGG_GSEA, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
View(head(KEGG_GSEA2))
x <- cnetplot(KEGG_GSEA2, categorySize="pvalue", foldChange=Book1_gsea, showCategory = KEGG_GSEA$Description[3])
x

x2 <- cnetplot(KEGG_GSEA2, categorySize="pvalue", foldChange=Book1_gsea, showCategory = 5, colorEdge = TRUE, layout = 'gem', cex_label_category = 0.5, cex_label_gene = 0.5)
x2 <- x2 + ggtitle("KEGG Category network")+
  theme(plot.title = element_text(size=15))
x2

tiff("./KEGG_cnetplot_5_cat.tiff", units = "in", height = 10, width = 10, res = 300)
x2
dev.off()

#Visualize first 5 GO pathways in cnetplot
GO_GSEA2 <- setReadable(GO_GSEA, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')
View(head(GO_GSEA2))

GO2 <- cnetplot(GO_GSEA2, categorySize="pvalue", foldChange=Book1_gsea, showCategory = 5, colorEdge = TRUE, layout = 'gem', cex_label_category = 0.8, cex_label_gene = 0.8)
GO2 <- GO2 + ggtitle("GO Category network")+
  theme(plot.title = element_text(size=15))
GO2

tiff("./GO_cnetplot_5_cat.tiff", units = "in", height = 15, width = 15, res = 300)
GO2
dev.off()
