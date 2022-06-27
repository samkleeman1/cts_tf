library(edgeR)
library(readr)
library(dplyr)
library(sva)

#Import raw reads (GRCh37-aligned, STAR counted)
raw <- read_csv("~/Dropbox (CSHL Dropbox Team)/PhD/Ben/CTS_RNA_seq/Oxford_CTS_raw_counts.csv")

genes = raw$...1
genes=sub(pattern="\\..*", replacement="", x=genes)

names = names(raw)[-c(1)]
batch = ifelse(grepl('sample',names)==TRUE,1,2)
raw=raw[,-c(1)]
adjusted <- ComBat_seq(as.matrix(raw), batch=batch, group=NULL) #Batch correction with ComBat-seq

#Add genome reference annotations
library(biomaRt)
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
#ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
conversion<-getBM(attributes = c("ensembl_gene_id", "external_gene_name","entrezgene_id"), 
                  filters = "ensembl_gene_id", 
                  values = genes, 
                  mart = ensembl)
n<-match(genes, conversion$ensembl_gene_id)
genes_full=data.frame(ensembl_id=genes, gene_name=conversion$external_gene_name[n], entrez_id = conversion$entrezgene_id[n])

library("GenomicFeatures")
gtf_txdb <- makeTxDbFromGFF("~/Dropbox (CSHL Dropbox Team)/Scratch/Downloads_old/gencode.v26.annotation.gtf.gz")
exons_list_per_gene <- exonsBy(gtf_txdb,by="gene")
widths <- width(reduce(exons_list_per_gene))
totalexonlength <- vapply(widths, sum, numeric(1))
m<-match(gsub("\\..*","",genes_full$ensembl_id), gsub("\\..*","",names(totalexonlength)))
genes_full$length = totalexonlength[m]

#Import counts data to edgeR
y <- DGEList(counts=as.matrix(adjusted),genes=genes_full)
keep <- rowSums(cpm(y)>1) >= 10 
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

#Extract log-CPM data (for inter-sample comparisons)
logcpm = cpm(y, log=TRUE)
row.names(logcpm) = y$genes$gene_name
logcpm = logcpm[!duplicated(row.names(logcpm)),]
logcpm = logcpm[!is.na(row.names(logcpm)),]
logcpm = logcpm[!is.na(rowSums(logcpm)),]

pca_res <- prcomp(t(logcpm), scale. = TRUE)
library(ggfortify)
autoplot(pca_res,label = TRUE) #Ensure no significant batch effect

#Pull in genomic data
logcpm=adj_logcpm
sub = logcpm[,grepl('sample', colnames(logcpm))==TRUE]
subx = logcpm[,grepl('sample', colnames(logcpm))==FALSE]

#Match RNA-seq IDs to genotyping IDs
sample_map = read_delim("~/Dropbox (CSHL Dropbox Team)/Scratch/Downloads_old/Translating sample to CTS_ID.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE)
sample_map$id=paste('sample',sample_map$wtchg_id, sep='')

m<-match(colnames(sub), sample_map$id)
colnames(sub)=sample_map$cts_id[m]

sub=cbind(sub, subx)

sub = subset(sub, row.names(sub)=="DIRC3" | row.names(sub)=="IGFBP5")
sub = as.data.frame(t(sub))

#Import genotype data
snp = read_csv("~/Dropbox (CSHL Dropbox Team)/Scratch/Downloads_old/cts_dirc3_snp.csv")
#Import principal component data
pcs = cts_pcs <- read_delim("~/Dropbox (CSHL Dropbox Team)/Scratch/Downloads_old/cts_pcs.tsv", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)

#Model data
m<-match(row.names(sub),snp$s)
sub$genotype = snp$score_allele[m]
m<-match(row.names(sub),pcs$IID)
sub$PC1 = pcs$PC1[m]
sub$PC2 = pcs$PC2[m]
sub$PC3 = pcs$PC3[m]
sub$PC4 = pcs$PC4[m]
m<-match(row.names(sub),sample_map$cts_id)
sub$sex=sample_map$gender[m]
sub$sex=ifelse(sub$sex==1,0,1)
sub = subset(sub, is.na(genotype)==FALSE)

model1 = lm(IGFBP5 ~ genotype + sex + PC1 + PC2 + PC3 + PC4, data=sub)

model2 = lm(DIRC3 ~ genotype + sex + PC1 + PC2 + PC3 + PC4, data=sub)


#Plot results
library(ggpubr)
p1 = ggboxplot(sub, 'genotype','DIRC3')+ggtitle('DIRC3')+ylim(-1,4)
p2 = ggboxplot(sub, 'genotype','IGFBP5')+ggtitle('IGFBP5')+ylim(7,13)
library(patchwork)
p1 + p2 + plot_layout(ncol = 2, heights = c(6)) + theme(panel.background = element_blank())



