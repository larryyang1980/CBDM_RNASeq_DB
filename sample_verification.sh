#!/bin/bash

#Purpose : Compare user selected sample with FASTQ directory to verify sample has corresponding FASTQ file
#Purpose : Create valid sample_list.txt file for use with pipeline.sh

#Output_1 : sample_list_report.txt listing samples with no corresponding FASTQ file
#Dependency : 'output' directory is setup to autosync with FileMaker Server/Data/Documents/ directory 
#Output_1 : sample_list_report.txt goes in output directory so it will sync back to FileMaker Documents in order to update database with
#Output_1 : non-matching samples

#Output_2 : sample_list.txt, a list of samples with corresponding FASTQ files that will be source of 
#Output_2 : pipeline.sh perform alignment script

#Absolute path to this script, e.g. /groups/cbdm-db/rnaseq_db/sample_verification.sh
#script=$(readlink -f "$0")
#Absolute path this script is in, e.g. /groups/cbdm-db/rnaseq_db
directory="/n/groups/cbdm-db/rnaseq_db"
#Get Run ID
RUN_ID=$1

#Clean up
rm -f $directory/Sample_List_Info/sample_list_report_${RUN_ID}.txt
rm -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt
rm -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt

#Make sure the file is in the right format
dos2unix $directory/Sample_List_Info/sample_list_unverified_${RUN_ID}.txt

#mkdir $directory/output/testoutput/

#Loop through the samples names

while IFS=$'\t' read -r -a Array
do
    samplename=${Array[0]}
    #facility sample name
    f=${Array[1]}
    owner=${Array[2]}
    species=${Array[3]}

#Loop through fastq files verifying sample name exists
for file in $directory/FASTQ/*; do
    #Filter on files
    if [[ -f $file ]]; then
      #Test if the file contains the sample name
      if [[ "$file" =~ "$f" ]]; then
         exist="true"
         echo -e $samplename'\t'$f'\t'$owner  >> $directory/Sample_List_Info/sample_list_${RUN_ID}_${species}.txt
         #echo -e $directory/output/$f >> $directory/output/sample_list_multiqc.txt
         break
      else
         exist="false"
      fi
    fi
done

#Creating the report for any  missing files
if [[ $exist == "false" ]]; then
   echo -e $samplename'\t'$f'\t'$owner'\t'$species  >> $directory/Sample_List_Info/sample_list_report_${RUN_ID}.txt
fi

done < $directory/Sample_List_Info/sample_list_unverified_${RUN_ID}.txt

echo "Launching the pipeline script..."
if [[ -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Homosapiens.txt ]]; then
	$directory/scripts/pipeline_O2.sh $directory/scripts/config_O2_human.txt $RUN_ID Homosapiens
fi

if [[ -f $directory/Sample_List_Info/sample_list_${RUN_ID}_Musmusculus.txt ]]; then
	$directory/scripts/pipeline_O2.sh $directory/scripts/config_O2_mouse.txt $RUN_ID Musmusculus
fi
