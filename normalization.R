if(!require(rmarkdown)){
	install.packages("rmarkdown", repos = "https://cloud.r-project.org/")
	library(rmarkdown)
}
if(!require(Hmisc)){
	install.packages("Hmisc", repos = "https://cloud.r-project.org/")
	library(Hmisc)
}
if(!require(DESeq2)){
	source("https://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
	library(DESeq2)
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
if(!require(matrixStats)){
	install.packages("matrixStats", repos = "https://cloud.r-project.org/")
	library(matrixStats)
}


args <- commandArgs(trailingOnly = TRUE)
manifest.file <- args[1]
RUN_ID <- args[2]
SPECIES <- args[3]

if ( length(args) < 3) {
    print ("Please input <sample info file> <RUN_ID> <SPECIES>\n")
    q(save = "no")
}

################################################################################

#norm_processing <- function(owner_id, owner_status, RUN_ID) {
#	multiqc_star_all <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output_tmp/normalization/",RUN_ID,"/multiqc_star.txt",sep=""),header=TRUE, comment.char = "", check.names = F, as.is = T)
#	multiqc_featureCounts_all <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output_tmp/normalization/",RUN_ID,"/multiqc_featureCounts.txt",sep=""),header=TRUE, comment.char = "", check.names = F, as.is = T)
#	if(owner_status) {
#		sample_info_all <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/sample_list.txt",sep=""), header=FALSE, stringsAsFactors = FALSE, comment.char = "", check.names = F, as.is = T)
#	} else {
#		sample_info_all <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/sample_list.txt",sep=""), header=FALSE, stringsAsFactors = FALSE, comment.char = "", check.names = F, as.is = T)
#	}
#	colnames(sample_info_all) <- c("Sample_Name", "Sample", "Owner")
#	rownames(sample_info_all) <- sample_info_all$Sample_Name
#	
#	if(owner_status) {
#		out_tbl <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/genes_count_table.tsv",sep=""),header=TRUE, comment.char = "", check.names = F, as.is = T)
#	} else {
#		out_tbl <- read.table(file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/genes_count_table.tsv",sep=""),header=TRUE, comment.char = "", check.names = F, as.is = T)
#	}
#	sample_names <-  data.frame(Sample_Name=colnames(out_tbl[,-1]), stringsAsFactors = FALSE)
#	
#	sample_info <- merge(sample_names, sample_info_all)
#	multiqc_star <- merge(sample_info, multiqc_star_all)
#	multiqc_featureCounts <- merge(sample_info, multiqc_featureCounts_all)
#	
#	multiqc_star_output <- multiqc_star[, c("Sample_Name","Sample","total_reads","uniquely_mapped","uniquely_mapped_percent")]
#	multiqc_featureCounts_output <- multiqc_featureCounts[, c("Sample_Name", "Sample", "Assigned", "percent_assigned", "Unassigned_Ambiguity", "Unassigned_MultiMapping", "Unassigned_NoFeatures")]
#	
#	sample_names_passing <- as.vector(multiqc_star_output$Sample_Name[multiqc_star_output$uniquely_mapped>1000000])
#	sample_names_all <- as.vector(multiqc_star_output$Sample_Name)
#	
#	counts <- out_tbl[, sample_names_passing]
#	counts_all <- out_tbl[, sample_names_all]
#	counts_output <- data.frame(out_tbl[,1], counts, stringsAsFactors = FALSE)
#	colnames(counts_output) = c("Gene_Symbol", colnames(counts))
#	if(owner_status) {
#		write.table(counts_output, file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/genes_count_table_QC_passed.tsv",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
#	} else {
#		write.table(counts_output, file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/genes_count_table_QC_passed.tsv",sep=""), quote=FALSE, row.names=FALSE, sep="\t")
#	}
#	
#	gene_names <- out_tbl[,1]
#	rownames(counts) <- gene_names
#	rownames(counts_all) <- gene_names
#	
#	num_of_genes_at_least_1_hit <- as.data.frame(sapply(out_tbl[,sample_names_all], function(x)sum(x>0)))
#	num_of_genes_over_10_hits <- as.data.frame(sapply(out_tbl[,sample_names_all], function(x)sum(x>10)))
#	num_of_genes_over_50_hits <- as.data.frame(sapply(out_tbl[,sample_names_all], function(x)sum(x>50)))
#	num_of_genes_mapped <- cbind(sample_names_all, num_of_genes_at_least_1_hit, num_of_genes_over_10_hits, num_of_genes_over_50_hits)
#	colnames(num_of_genes_mapped) <- c("Sample_Name", "Number of genes with at least 1 mapped reads", "Number of genes with over 10 mapped reads", "Number of genes with over 50 mapped reads")
#	
#	if(owner_status) {
#		write.xlsx2(merge(merge(multiqc_star_output,multiqc_featureCounts_output),num_of_genes_mapped), file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/Reads_Mapping_Quantification_Statistics.xlsx",sep=""), sheetName = "Statistics", row.names = F)
#	} else {
#		write.xlsx2(merge(merge(multiqc_star_output,multiqc_featureCounts_output),num_of_genes_mapped), file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/Reads_Mapping_Quantification_Statistics.xlsx",sep=""), sheetName = "Statistics", row.names = F)
#	}
#	
#	
#	Num_of_samples <- length(sample_names_passing)
#	sample_of_class <- unique(unlist(lapply(sample_names_passing, function(x) strsplit(x, split="#")[[1]][1])))
#	Num_of_class <- length(sample_of_class)
#	sample_index <- matrix(match(unlist(lapply(sample_names_passing, function(x) strsplit(x, split="#")[[1]][1])), sample_of_class)-1,nrow=1,byrow=TRUE)
#	cls_header_1 <- matrix(c(Num_of_samples, Num_of_class, 1),nrow=1,byrow=TRUE)
#	cls_header_2 <- matrix(c('#', sample_of_class),nrow=1,byrow=TRUE)
#	if(owner_status) {
#		write.table(cls_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE)
#		write.table(cls_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE,append=TRUE)
#		write.table(sample_index, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.cls",sep=""),sep=' ',quote=FALSE, row.names=FALSE,col.names=FALSE,append=TRUE)
#	} else {
#		write.table(cls_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE)
#		write.table(cls_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.cls",sep=""),quote=FALSE,sep = ' ',row.names=FALSE,col.names=FALSE,append=TRUE)
#		write.table(sample_index, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.cls",sep=""),sep=' ',quote=FALSE, row.names=FALSE,col.names=FALSE,append=TRUE)
#	}
#	
#	pData <- sample_info_all[sample_names_passing,]
#	deseqobj <- DESeqDataSetFromMatrix(countData=counts, colData=pData, design= ~ 1)
#	dds <- estimateSizeFactors(deseqobj)
#	normalized.counts <-counts(dds, normalized=TRUE)
#	
#	dat <- normalized.counts
#	dat <- dat+1
#	pt1_adjust <- data.frame(gene_names,gene_names,dat)
#	col_names <- c('NAME','DESCRIPTION',sample_names_passing)
#	names(pt1_adjust) <- col_names
#	#write out a .gct using the appropriate format
#	gct_header_1 <- matrix(c('#1.2'),nrow=1,byrow=TRUE)
#	gct_header_2 <- matrix(c(dim(dat)[1],dim(dat)[2]),nrow=1,byrow=TRUE)
#	if(owner_status) {
#		write.table(gct_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE)
#		write.table(gct_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE,append=TRUE)
#		write.table(pt1_adjust,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_1_adjust.gct",sep=""),sep='\t',quote=FALSE,row.names=FALSE,append=TRUE)
#	} else {
#		write.table(gct_header_1, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE)
#		write.table(gct_header_2, file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.gct",sep=""),quote=FALSE,sep = '\t',row.names=FALSE,col.names=FALSE,append=TRUE)
#		write.table(pt1_adjust,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_1_adjust.gct",sep=""),sep='\t',quote=FALSE,row.names=FALSE,append=TRUE)
#	}
#	
#	pData_all <- sample_info_all[sample_names_all,]
#	deseqobj_all <- DESeqDataSetFromMatrix(countData=counts_all, colData=pData_all, design= ~ 1)
#	dds_all <- estimateSizeFactors(deseqobj_all)
#	normalized.counts_all <-counts(dds_all, normalized=TRUE)
#	
#	#calculate mean and sd of expression across genes
#	gene_means <- rowMeans(normalized.counts_all)
#	gene_sds <- apply(normalized.counts_all,1,sd)
#	gene_cv <- gene_sds/gene_means
#	#calculate mean and sd of expression across samples
#	sample_means <- colMeans(normalized.counts_all)
#	sample_sds <- apply(normalized.counts_all,2,sd)
#	sample_cv <- sample_sds/sample_means
#	#calculate gene expression data
#	gene_expression_data <- data.frame(gene_names,gene_means,gene_sds,gene_cv)
#	names(gene_expression_data) <- c('gene_names','gene_means','gene_sds','gene_cv')
#	
#	sample_data <- data.frame(names(normalized.counts_all[1,]),sample_means,sample_sds,sample_cv)
#	names(sample_data) <- c('sample_names','sample_means','sample_sds','sample_cv')
#	write.table(sample_data,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/sample_expression_data.csv",sep=""),quote=FALSE,row.names=FALSE,sep=',')
#	
#	#calculate Pearson correlation matrix across samples
#	corr_data <- cor(normalized.counts_all)
#	corr_data_output <- data.frame(rownames(corr_data),corr_data)
#	colnames(corr_data_output) = c("Sample_Names", colnames(corr_data))
#	
#	if(owner_status) {
#		write.table(gene_expression_data,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/gene_expression_data.csv",sep=""),row.names=FALSE,quote=FALSE,sep=',')
#		write.table(sample_data,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/sample_expression_data.csv",sep=""),quote=FALSE,row.names=FALSE,sep=',')
#		write.table(corr_data_output, quote=FALSE, row.names=FALSE, file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,"/Pearson_Correlation_Matrix_Samples.csv",sep=""), sep=',')
#	} else {
#		write.table(gene_expression_data,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/gene_expression_data.csv",sep=""),row.names=FALSE,quote=FALSE,sep=',')
#		write.table(sample_data,file=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/sample_expression_data.csv",sep=""),quote=FALSE,row.names=FALSE,sep=',')
#		write.table(corr_data_output, quote=FALSE, row.names=FALSE, file = paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/Pearson_Correlation_Matrix_Samples.csv",sep=""), sep=',')
#	}
#	
#	if(owner_status) {
#		rmarkdown::render("/n/groups/cbdm-db/rnaseq_db/RNASeq_QC_Report.Rmd", output_dir=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",owner_id,sep=""), params=list(tag=owner_id))
#	} else {
#		rmarkdown::render("/n/groups/cbdm-db/rnaseq_db/RNASeq_QC_Report.Rmd", output_dir=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,sep=""), params=list(tag=RUN_ID))
#	}
#}

################################################################################

#read all sample id
manifest <- read.table(manifest.file, header=F, sep="\t", comment.char="", stringsAsFactors=F)

Sample_Names <- manifest[,1]
Sample_IDs <- manifest[,2]
Owners <- manifest[,3]

results = list()
i=0
for (id in Sample_IDs) {
    filename = paste("/n/groups/cbdm-db/rnaseq_db/output/alignment/", RUN_ID, "/", SPECIES, "/Count_Tables/", id, "/", id, ".featureCounts.txt", sep="")

    data <- read.table(filename, header=TRUE, stringsAsFactors=F, check.names = F)
    data <- data[with(data, order(Geneid)),]

    results[[id]] = data[,7]
    if ( i < 1 ) {
	    genes=data[,1]
	    gene_lens=data[,6]
    }
    i <- i+1
}

counts= as.data.frame( do.call("cbind", results))
expr.count = cbind( genes, counts )
colnames(expr.count) <- c("Gene_Symbol",Sample_Names)

# write raw read counts table
write.table(expr.count, paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",SPECIES,"/genes_count_table.tsv",sep=""), sep="\t",row.names = F, quote = F)
#norm_processing(RUN_ID,0,RUN_ID)
rmarkdown::render("/n/groups/cbdm-db/rnaseq_db/scripts/RNASeq_QC_Report.Rmd", output_dir=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",SPECIES,sep=""), params=list(owner_id=RUN_ID,owner_status=0,RUN_ID=RUN_ID,SPECIES=SPECIES))

for(owner in unique(Owners)) {
	samples_tmp = Sample_Names[Owners == owner]
	expr.count_tmp = expr.count[,c("Gene_Symbol",samples_tmp)]
	write.table(expr.count_tmp, paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",SPECIES,"/",owner,"/genes_count_table.tsv",sep=""), sep="\t",row.names = F, quote = F)
	#norm_processing(owner,1,RUN_ID)
	rmarkdown::render("/n/groups/cbdm-db/rnaseq_db/scripts/RNASeq_QC_Report.Rmd", output_dir=paste("/n/groups/cbdm-db/rnaseq_db/output/normalization/",RUN_ID,"/",SPECIES,"/",owner,sep=""), params=list(owner_id=owner,owner_status=1,RUN_ID=RUN_ID,SPECIES=SPECIES))

}
