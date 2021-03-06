#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

DOCKER=$( getparameter docker )
NAME=$( getparameter name )
GENOME=$( getparameter genome )

if [ -n "$DOCKER" ]; then
    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
fi

set -ex

if [ -n "${input_bed}" ]; then
    BED=$(tempfile)
    tail -n +2 ${input_bed} | sed 's/"//g' > $BED
else
	if [ "$GENOME" = "human/hg19" ]; then
		BED=/home/anduril/BED/Human_Homo_sapiens/hg19.HouseKeepingGenes.nochr.bed
	else
		BED=/home/anduril/BED/Mouse_Mus_musculus/mm10_Ensembl.nochr.bed
	fi
fi

cd $(dirname ${output_coverage})
$DOCKER geneBody_coverage.py \
	-i ${input_bam} \
	-r $BED \
	-o $(dirname ${output_coverage})/output

cat $(dirname ${output_coverage})/output.geneBodyCoverage.txt | perl -ne "if (!/Percentile/) { s/^\S+/$NAME/ }; print \$_" > ${output_coverage}	
rm $(dirname ${output_coverage})/output.geneBodyCoverage.txt
