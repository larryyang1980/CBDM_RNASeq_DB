if(!require(rmarkdown)){
	install.packages("rmarkdown", repos = "https://cloud.r-project.org/")
	library(rmarkdown)
}
if(!require(Hmisc)){
	install.packages("Hmisc", repos = "https://cloud.r-project.org/")
	library(Hmisc)
}
if(!require(dplyr)){
	install.packages("dplyr", repos = "https://cloud.r-project.org/")
	library(dplyr)
}
if(!require(readr)){
	install.packages("readr", repos = "https://cloud.r-project.org/")
	library(readr)
}
if(!require(knitr)){
	install.packages("knitr", repos = "https://cloud.r-project.org/")
	library(knitr)
}
if(!require(DT)){
	install.packages("DT", repos = "https://cloud.r-project.org/")
	library(DT)
}
if(!require(heatmaply)){
	install.packages("heatmaply", repos = "https://cloud.r-project.org/")
	library(heatmaply)
}
if(!require(plotly)){
	install.packages("plotly", repos = "https://cloud.r-project.org/")
	library(plotly)
}
if(!require(xlsx)){
	install.packages("xlsx", repos = "https://cloud.r-project.org/")
	library(xlsx)
}
if(!require(DESeq2)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
	library(DESeq2)
}

args <- commandArgs(trailingOnly = TRUE)
manifest.file <- args[1]
USER_ID <- args[2]
TIMESTAMP <- args[3]

if ( length(args) < 3) {
    print ("Please input <sample info file> <USER_ID> <TIMESTAMP>\n")
    q(save = "no")
}

################################################################################
## Generating Gene Counts Table
################################################################################

#read all sample id
manifest <- read.table(manifest.file, header=F, sep="\t", comment.char="", stringsAsFactors=F, check.names = F, as.is = T)

Sample_Names <- manifest[,1]
Sample_IDs <- manifest[,2]
RUN_IDs <- manifest[,3]
SPECIES <- manifest[,4]

results = list()
for (i in 1:length(Sample_IDs)) {
    filename = paste("/n/groups/cbdm-db/rnaseq_db/output/alignment/", RUN_IDs[i], "/", SPECIES[i], "/Count_Tables/", Sample_IDs[i], "/", Sample_IDs[i], ".featureCounts.txt", sep="")

    data <- read.table(filename, header=TRUE, stringsAsFactors=F, check.names = F)
    data <- data[with(data, order(Geneid)),]

    results[[Sample_IDs[i]]] = data[,7]
    if ( i == 1 ) {
	    genes=data[,1]
	    #gene_lens=data[,6]
    }
}

counts= as.data.frame( do.call("cbind", results))
expr.count = cbind( genes, counts )
colnames(expr.count) <- c("Gene_Symbol",Sample_Names)

# write raw read counts table
write.table(expr.count, paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/genes_count_table_",TIMESTAMP,".tsv",sep=""), sep="\t",row.names = F, quote = F)


################################################################################
## Normalization
################################################################################
sample_info <- read.table(manifest.file, header=F, sep="\t", comment.char="", stringsAsFactors=F, check.names = F, as.is = T)
colnames(sample_info) <- c("Sample_Name", "Sample", "Batch", "Species")
rownames(sample_info) <- sample_info$Sample_Name

out_tbl <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/genes_count_table_",TIMESTAMP,".tsv",sep=""),header=TRUE, comment.char = "", check.names = F, as.is = T)
counts <- out_tbl[,-1]
gene_names <- out_tbl[,1]
rownames(counts) <- gene_names

sample_names_passing = sample_info$Sample_Name

Num_of_samples <- length(sample_names_passing)
sample_of_class <- unique(unlist(lapply(sample_names_passing, function(x) strsplit(x, split="#")[[1]][1])))
Num_of_class <- length(sample_of_class)
sample_index <- matrix(match(unlist(lapply(sample_names_passing, function(x) strsplit(x, split="#")[[1]][1])), sample_of_class)-1,nrow=1,byrow=TRUE)
cls_header_1 <- matrix(c(Num_of_samples, Num_of_class, 1),nrow=1,byrow=TRUE)
cls_header_2 <- matrix(c('#', sample_of_class),nrow=1,byrow=TRUE)
write.table(cls_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE)
write.table(cls_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(sample_index, file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".cls",sep=""),sep=' ',quote=FALSE, row.names=FALSE,col.names=FALSE,append=TRUE)

pData <- sample_info[sample_names_passing,]
if(nlevels(pData$Batch) == 1) {
	deseqobj <- DESeqDataSetFromMatrix(countData=counts, colData=pData, design= ~ 1)
} else {
	deseqobj <- DESeqDataSetFromMatrix(countData=counts, colData=pData, design= ~ Batch)
}
dds <- estimateSizeFactors(deseqobj)
normalized.counts <-counts(dds, normalized=TRUE)

data_1_adjust <- data.frame(gene_names,gene_names,normalized.counts+1)
names(data_1_adjust) <- c('NAME','DESCRIPTION',sample_names_passing)
gct_header_1 <- matrix(c('#1.2'),nrow=1,byrow=TRUE)
gct_header_2 <- matrix(c(dim(normalized.counts)[1],dim(normalized.counts)[2]),nrow=1,byrow=TRUE)
write.table(gct_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE)
write.table(gct_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(data_1_adjust,file=paste("/n/groups/cbdm-db/rnaseq_db/output/analyze_datasets/",USER_ID,"/normalized_gene_expression_1_adjust_",TIMESTAMP,".gct",sep=""),sep='\t',quote=FALSE,row.names=FALSE,append=TRUE)
