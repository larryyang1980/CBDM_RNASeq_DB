#!/bin/bash

module load gcc/6.2.0
module load R/3.4.1

#Absolute path to this script, e.g. /groups/cbdm-db/rnaseq_db/normalization.sh
#script=$(readlink -f "$0")
#Absolute path this script is in, e.g. /groups/cbdm-db/rnaseq_db
directory="/n/groups/cbdm-db/rnaseq_db"
SCRIPT_DIR="/n/groups/cbdm-db/ly82/RNASeq_pipeline/STAR_index"
LOG_DIR="/n/groups/cbdm-db/rnaseq_db/logs"

RUN_ID=$1

mkdir -p $directory/output_tmp/normalization
rm -fr $directory/output_tmp/normalization/$RUN_ID
mkdir -p $directory/output_tmp/normalization/$RUN_ID
mkdir -p $directory/output/normalization
rm -fr $directory/output/normalization/$RUN_ID
mkdir -p $directory/output/normalization/$RUN_ID
cd $directory/output/normalization/$RUN_ID

if [[ -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt ]]; then
	mkdir -p $directory/output/normalization/$RUN_ID/Homosapiens
	mkdir -p $directory/output_tmp/normalization/$RUN_ID/Homosapiens
	cp $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt $directory/output/normalization/$RUN_ID/Homosapiens/sample_list.txt
	cp $directory/output/alignment/$RUN_ID/Homosapiens/QC_data/*.txt $directory/output_tmp/normalization/$RUN_ID/Homosapiens/
	
	while IFS=$'\t' read -r -a Array
	do
	    Sample_Name=${Array[0]}
	    Sample_ID=${Array[1]}
	    Owner=${Array[2]}
	
	    mkdir -p $directory/output/normalization/$RUN_ID/Homosapiens/$Owner
	    
	    echo -e $Sample_Name'\t'$Sample_ID'\t'$Owner  >> $directory/output/normalization/$RUN_ID/Homosapiens/$Owner/sample_list.txt
	
	done < $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt
	
	
	###############################################################################
	### Generating Genes Count Table
	###############################################################################
	rm -f $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	cp $SCRIPT_DIR/multiqc.sh $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	cd $directory/output_tmp/normalization/$RUN_ID/Homosapiens
	
	echo -e "#SBATCH -c 4" >> $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	echo -e "#SBATCH -o $LOG_DIR/${RUN_ID}_Homosapiens.sbatch_normalization.log" >> $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	echo -e "module load gcc/6.2.0" >> $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	echo -e "module load R/3.4.1" >> $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	echo -e "Rscript $directory/scripts/normalization.R $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt $RUN_ID Homosapiens" >> $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
	
	sbatch $directory/output_tmp/normalization/$RUN_ID/Homosapiens/sbatch_normalization_Homosapiens.sh
fi

if [[ -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt ]]; then
	mkdir -p $directory/output/normalization/$RUN_ID/Musmusculus
	mkdir -p $directory/output_tmp/normalization/$RUN_ID/Musmusculus
	cp $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt $directory/output/normalization/$RUN_ID/Musmusculus/sample_list.txt
	cp $directory/output/alignment/$RUN_ID/Musmusculus/QC_data/*.txt $directory/output_tmp/normalization/$RUN_ID/Musmusculus/
	
	while IFS=$'\t' read -r -a Array
	do
	    Sample_Name=${Array[0]}
	    Sample_ID=${Array[1]}
	    Owner=${Array[2]}
	
	    mkdir -p $directory/output/normalization/$RUN_ID/Musmusculus/$Owner
	    
	    echo -e $Sample_Name'\t'$Sample_ID'\t'$Owner  >> $directory/output/normalization/$RUN_ID/Musmusculus/$Owner/sample_list.txt
	
	done < $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt
	
	
	###############################################################################
	### Generating Genes Count Table
	###############################################################################
	rm -f $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	cp $SCRIPT_DIR/multiqc.sh $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	cd $directory/output_tmp/normalization/$RUN_ID/Musmusculus
	
	echo -e "#SBATCH -c 4" >> $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	echo -e "#SBATCH -o $LOG_DIR/${RUN_ID}_Musmusculus.sbatch_normalization.log" >> $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	echo -e "module load gcc/6.2.0" >> $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	echo -e "module load R/3.4.1" >> $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	echo -e "Rscript $directory/scripts/normalization.R $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt $RUN_ID Musmusculus" >> $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
	
	sbatch $directory/output_tmp/normalization/$RUN_ID/Musmusculus/sbatch_normalization_Musmusculus.sh
fi


