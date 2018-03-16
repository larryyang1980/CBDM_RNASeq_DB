#!/bin/bash

module load gcc/6.2.0
module load R/3.4.1

#Absolute path to this script, e.g. /groups/cbdm-db/rnaseq_db/analyze_datasets.sh
#script=$(readlink -f "$0")
#Absolute path this script is in, e.g. /groups/cbdm-db/rnaseq_db
directory="/n/groups/cbdm-db/rnaseq_db"
SCRIPT_DIR="/n/groups/cbdm-db/ly82/RNASeq_pipeline/STAR_index"
LOG_DIR="/n/groups/cbdm-db/rnaseq_db/logs"

USER_ID=$1
TIMESTAMP=$2

mkdir -p $directory/output_tmp/analyze_datasets
mkdir -p $directory/output_tmp/analyze_datasets/$USER_ID
mkdir -p $directory/output/analyze_datasets
mkdir -p $directory/output/analyze_datasets/$USER_ID
cd $directory/output/analyze_datasets/$USER_ID
cp $directory/Sample_List_Info/sample_list_${USER_ID}_${TIMESTAMP}.txt $directory/output/analyze_datasets/$USER_ID/

###############################################################################
### Generating Genes Count Table
###############################################################################
rm -f $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
cp $SCRIPT_DIR/multiqc.sh $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
cd $directory/output_tmp/analyze_datasets/$USER_ID

echo -e "#SBATCH -c 4" >> $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
echo -e "#SBATCH -o $LOG_DIR/${USER_ID}.sbatch_analyze_datasets.log" >> $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
echo -e "module load gcc/6.2.0" >> $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
echo -e "module load R/3.4.1" >> $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh
echo -e "Rscript $directory/scripts/analyze_datasets.R $directory/Sample_List_Info/sample_list_${USER_ID}_${TIMESTAMP}.txt $USER_ID $TIMESTAMP" >> $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh

sbatch $directory/output_tmp/analyze_datasets/$USER_ID/sbatch_analyze_datasets.sh

