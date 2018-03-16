############################################################################
##
##  1. reads mapping
##  2. quantify gene expression level
##  3. aligned data QC
##  
##  Usage:
##  	pipeline.sh configure-file RUN_ID
##
##  Liang Yang (lyang@hms.harvard.edu)
##
##  February 27, 2017
##  Last modified:  February 27, 2017
##
############################################################################

#!/bin/bash

if [ $# -ne 2 ]; then
        echo "Usage: $0 <Configure file> <RUN_ID>"
	exit
fi

script=$(readlink -f "$0")
directory=$(dirname "$script")

CONFIG=$1
RUN_ID=$2

if [ -f $CONFIG ]
then
	#source $CONFIG
	#FASTQ_DIR=$(readlink -f $FASTQ_DIR)
	### Remove tailing white space...it screws up the naming schemes
	sed -e 's/ *$//' $CONFIG > tmp.txt 
	mv tmp.txt $CONFIG
	dos2unix $CONFIG

	#FLIST=`grep -w '^FLIST' $CONFIG | cut -d '=' -f2`
	FASTQ_DIR=`grep -w '^FASTQ_DIR' $CONFIG | cut -d '=' -f2`
	FASTQ_SUFFIX=`grep -w '^FASTQ_SUFFIX' $CONFIG | cut -d '=' -f2`
	SEQUENCE_TYPE=`grep -w '^SEQUENCE_TYPE' $CONFIG | cut -d '=' -f2`
	STRAND=`grep -w '^STRAND' $CONFIG | cut -d '=' -f2`
	Mapping_software=`grep -w '^Mapping_software' $CONFIG | cut -d '=' -f2`
	Count_software=`grep -w '^Count_software' $CONFIG | cut -d '=' -f2`
	QC_before_mapping_software=`grep -w '^QC_before_mapping_software' $CONFIG | cut -d '=' -f2`
	QC_software=`grep -w '^MultiQC' $CONFIG | cut -d '=' -f2`
	MultiQC=`grep -w '^QC_software' $CONFIG | cut -d '=' -f2`
	Samtools=`grep -w '^Samtools' $CONFIG | cut -d '=' -f2`
	LOG_DIR=`grep -w '^LOG_DIR' $CONFIG | cut -d '=' -f2`
	OUTPUT_DIR=`grep -w '^OUTPUT_DIR' $CONFIG | cut -d '=' -f2`
	OUTPUT_TMP_DIR=`grep -w '^OUTPUT_TMP_DIR' $CONFIG | cut -d '=' -f2`
	FASTQC_DIR=`grep -w '^FASTQC_DIR' $CONFIG | cut -d '=' -f2`
	GENOME_INDEX=`grep -w '^GENOME_INDEX' $CONFIG | cut -d '=' -f2`
	GTF_FILE=`grep -w '^GTF_FILE' $CONFIG | cut -d '=' -f2`
	BED_FILE=`grep -w '^BED_FILE' $CONFIG | cut -d '=' -f2`
	Mapping_PARAMETER=`grep -w '^Mapping_PARAMETER' $CONFIG | cut -d '=' -f2`
	FEATURECOUNTS_PARAMETER=`grep -w '^FEATURECOUNTS_PARAMETER' $CONFIG | cut -d '=' -f2`
	
else
	echo "Please provide a configuration file"
	exit
fi

###############################################################################
### Set Enviroments
###############################################################################
module load $Mapping_software
module load $Count_software
module load $QC_before_mapping_software
module load $QC_software
module load $Samtools
module load $MultiQC

#rm -fr $OUTPUT_TMP_DIR
mkdir -p $LOG_DIR
mkdir -p $OUTPUT_DIR
mkdir -p $OUTPUT_DIR/$RUN_ID
mkdir -p $OUTPUT_DIR/$RUN_ID/BAMs
mkdir -p $OUTPUT_DIR/$RUN_ID/Count_Tables
mkdir -p $OUTPUT_DIR/$RUN_ID/QC_data
mkdir -p $OUTPUT_TMP_DIR
#mkdir -p $FASTQC_DIR

# fq.gz or fastq.gz or fastq  or fq
READCMD="--readFilesCommand zcat"
if [[ $FASTQ_SUFFIX = "fastq" || $FASTQ_SUFFIX = "fq" ]]; then
     READCMD=""
fi

#To handle very deep sequencing samples
CORE=4
#limitBAMsortRAM=""
#MEMORY_USAGE=""
#M_PARAMETER=""
#if [[ $SEQUENCE_DEPTH = "deep" ]]; then
#	CORE=15
#        MEMORY_USAGE="rusage [mem=125829120]"
#	M_PARAMETER="-M 125829120"
#	limitBAMsortRAM="--limitBAMsortRAM 107374182400"
#fi

###############################################################################
### Mapping
###############################################################################
while IFS=$'\t' read -r -a Array
do
   f=${Array[1]}
   echo -e $OUTPUT_TMP_DIR/$f >> $OUTPUT_TMP_DIR/sample_list_multiqc_${RUN_ID}.txt
   
   mkdir -p $OUTPUT_TMP_DIR/$f
   mkdir -p $OUTPUT_DIR/$RUN_ID/BAMs/$f
   mkdir -p $OUTPUT_DIR/$RUN_ID/Count_Tables/$f
   cd $OUTPUT_TMP_DIR/$f
   
   # handle pair or single RNA-seq
	FASTQS="$FASTQ_DIR/${f}.R1.$FASTQ_SUFFIX $FASTQ_DIR/${f}.R2.$FASTQ_SUFFIX"
	PAIRENDS="-p -B -C"
	if [[ $SEQUENCE_TYPE = "single" ]]; then
		FASTQS="$FASTQ_DIR/${f}.$FASTQ_SUFFIX"
		PAIRENDS=""
	fi
	
	#bsub -q short -W 12:00 -J ${f}.fastqc -o $LOG_DIR/${f}.fastqc.log "fastqc -o $FASTQC_DIR $FASTQS"
	bsub -q short -W 12:00 -J ${f}.mapping -R "rusage[mem=50000]" -o $LOG_DIR/${f}.mapping.log "STAR --genomeDir $GENOME_INDEX --sjdbGTFfile $GTF_FILE --readFilesIn $FASTQS $READCMD --runThreadN $CORE --outSAMtype BAM SortedByCoordinate $Mapping_PARAMETER; mv Aligned.sortedByCoord.out.bam $f.sortedByCoord.bam; samtools index $f.sortedByCoord.bam; rm -rf _STARtmp _STARgenome"
	bsub -q short -W 12:00 -J ${f}.count -w "ended(${f}.mapping)" -o $LOG_DIR/${f}.count.log "featureCounts $PAIRENDS -T 4 -F GTF -a $GTF_FILE -g gene_name -s $STRAND $FEATURECOUNTS_PARAMETER -o ${f}.featureCounts.txt $f.sortedByCoord.bam"

done < $directory/Sample_List_Info/sample_list_${RUN_ID}.txt

###############################################################################
### MultiQC
###############################################################################
cd $OUTPUT_TMP_DIR
bsub -q short -W 12:00 -J ${RUN_ID}.multiqc -w "ended(*.count)" -o $LOG_DIR/${RUN_ID}.multiqc.log "multiqc -f -l -n multiqc_${RUN_ID} $OUTPUT_TMP_DIR/sample_list_multiqc_${RUN_ID}.txt"

###############################################################################
### Copying files to output folder
###############################################################################
bsub -q short -W 12:00 -J ${RUN_ID}.copying -w "ended(${RUN_ID}.multiqc)" -o $LOG_DIR/${RUN_ID}.copying.log "/groups/cbdm-db/rnaseq_db/copy_files.sh $directory/Sample_List_Info/sample_list_${RUN_ID}.txt $OUTPUT_DIR $OUTPUT_TMP_DIR $RUN_ID"

