#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

DOCKER=$( getparameter docker )
NAME=$( getparameter name )

if [ -n "$DOCKER" ]; then
    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
fi

set -ex

$DOCKER geneBody_coverage.py \
	-i ${input_bam} \
	-r /home/anduril/BED/Human_Homo_sapiens/hg19.HouseKeepingGenes.nochr.bed \
	-o $(dirname ${output_coverage})/output

cat $(dirname ${output_coverage})/output.geneBodyCoverage.txt | perl -ne "if (!/Percentile/) { s/^\S+/$NAME/ }; print \$_" > ${output_coverage}	
rm $(dirname ${output_coverage})/output.geneBodyCoverage.txt
