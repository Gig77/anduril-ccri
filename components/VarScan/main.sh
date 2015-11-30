#!/bin/bash

source "$ANDURIL_HOME/bash/functions.sh"
export_command

set -ex
set -o pipefail

java -XX:+UseParallelGC -XX:ParallelGCThreads=8 -jar ${parameter_jar} mpileup2cns \
	<(samtools mpileup ${parameter_optionsSam} -l <(zcat /mnt/projects/iamp/data/Homo_sapiens.GRCh37.75.etv6runx1.gtf.gz | grep -P "\tCDS\t" | cut -f 1,4,5) -f ${input_reference} ${input_bam1}) \
	${parameter_options} \
	> ${output_variants}.part

mv ${output_variants}.part ${output_variants}

