#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

set -ex

#cut -f 2 ${input__index_bamFiles} | tail -n +2 > /root/bamlist.txt

geneBody_coverage.py \
	-i ${input_bam} \
	-r /root/BED/Human_Homo_sapiens/hg19.HouseKeepingGenes.nochr.bed \
	-o $(dirname ${output_coverage})/output
	
mv $(dirname ${output_coverage})/output.geneBodyCoverage.txt ${output_coverage} 
