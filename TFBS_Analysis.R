if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("bedr")

setwd("~/Documents/Opleiding/Biomedical Sciences/master/Jaar 1/From Big data to personalised medicine/report/")

getwd()

DE <- read.table("sorted_DE_stim2_4.csv", header=T, sep=",")

rownames(DE) <- DE$X
DE <- DE[,-c(1)]

library(biomaRt)

ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)

genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(DE),
                 mart = ensembl )

DE_extended <- merge(x = DE, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(DE_extended)
prom=matrix(ncol=3, nrow=nrow(DE_extended))
for (i in 1:nrow(DE_extended)){
  if (DE_extended$strand[i] == "1") 
  {prom[i,1] <- as.character(DE_extended$Row.names[i])
  prom[i,2] <- DE_extended$start_position[i]-350
  prom[i,3] <- DE_extended$start_position[i]+150
  }
  else {
    prom[i,1] <- as.character(DE_extended$Row.names[i])
    prom[i,2] <- DE_extended$end_position[i]-150
    prom[i,3] <- DE_extended$end_position[i]+350}   
}
colnames(prom)=c("Row.names","prom_start","prom_end")
comb <- merge(DE_extended, prom, by.x="Row.names", by.y="Row.names")
DE250 <- comb[,c("chromosome_name","prom_start", "prom_end","Row.names","hgnc_symbol","strand")]
colnames(DE250) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
DE250$chr <- paste("chr", DE250$chr, sep="")
DE250$chr <- as.character(as.factor(DE250$chr))
DE250$start <- as.numeric(as.character(DE250$start))
DE250$end <- as.numeric(as.character(DE250$end))
DE250$ENSG <- as.character(DE250$ENSG)
CTR <- read.table("control_Genes.csv", header=T, sep=",", row.names = 2)

library(biomaRt)
ensembl  <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol","gene_biotype","chromosome_name","start_position","end_position","strand"),
                 filters = "ensembl_gene_id",
                 values = rownames(CTR),
                 mart = ensembl )
CTR_extended <- merge(x = CTR, y = genemap, by.x = 0, by.y = "ensembl_gene_id", all.x = TRUE, all.y = FALSE)
head(CTR_extended)


prom=matrix(ncol=3, nrow=nrow(CTR_extended))
for (i in 1:nrow(CTR_extended)){
  if (CTR_extended$strand[i] == "1") 
  {prom[i,1] <- as.character(CTR_extended$Row.names[i])
  prom[i,2] <- CTR_extended$start_position[i]-350
  prom[i,3] <- CTR_extended$start_position[i]+150
  }
  else {
    prom[i,1] <- as.character(CTR_extended$Row.names[i])
    prom[i,2] <- CTR_extended$end_position[i]-150
    prom[i,3] <- CTR_extended$end_position[i]+350}  
}
colnames(prom)=c("Row.names","prom_start","prom_end")
comb <- merge(CTR_extended, prom, by.x="Row.names", by.y="Row.names")
CTR250 <- comb[,c("chromosome_name","prom_start", "prom_end","Row.names","hgnc_symbol","strand")]
colnames(CTR250) <- c("chr","start", "end", "ENSG", "hgnc", "strand")
CTR250$chr <- paste("chr", CTR250$chr, sep="")
CTR250$chr <- as.character(as.factor(CTR250$chr))
CTR250$start <- as.numeric(as.character(CTR250$start))
CTR250$end <- as.numeric(as.character(CTR250$end))

library(bedr)
tmp <- as.data.frame(CTR250[,c(1:3)])
View(tmp)
CTRL_SORT <- bedr.sort.region(tmp)

tmp2 <- as.data.frame(DE250[,c(1:3)])
View(tmp2)
DE_SORT <- bedr.sort.region(tmp2)

DE_SORT <- merge(x= DE_SORT,y= DE250[,c('start','ENSG', 'hgnc', 'strand')], by = 'start', sort = FALSE)
DE_SORT = DE_SORT[,c(2,1,3,4,5,6)]

CTRL_SORT <- merge(x= CTRL_SORT,y= CTR250[,c('start','ENSG', 'hgnc', 'strand')], by = 'start', sort = FALSE)
CTRL_SORT = CTRL_SORT[,c(2,1,3,4,5,6)]

write.csv(DE_SORT, file = "~/Documents/Opleiding/Biomedical Sciences/master/Jaar 1/From Big data to personalised medicine/report/DE_Sorted.csv", row.names = FALSE) 

write.csv(CTRL_SORT, file = "~/Documents/Opleiding/Biomedical Sciences/master/Jaar 1/From Big data to personalised medicine/report/CTRL_Sorted.csv", row.names = FALSE) 
