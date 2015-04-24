#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"

export_command

set -ex

mkdir "${output_log}"

cd $( dirname ${output_alignment} )

# Copy parameter file to cwd so HISAT will use it
#if [[ "${parameter_parameters}" != "" ]]; then
#    echo "Using parameter_parameters from: ${parameter_parameters}"
#    cp ${parameter_parameters} Parameters1.in
#fi

READS=$( getinput reads )
echo "Reads are ${READS}"

OPTIONS=$( getparameter options )

# align
/home/anduril/hisat-0.1.5-beta/hisat \
	-x /home/anduril/index/human_g1k_v37 \
	-U <(java -jar /home/anduril/picard-tools-1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${READS} FASTQ=/dev/stdout) \
	${OPTIONS} \
		| /home/anduril/samtools-1.2/samtools view -Shb - \
		> ${output_alignment}.unsorted

# sort
/home/anduril/samtools-1.2/samtools sort ${output_alignment}.unsorted ${output_alignment}.sorted
rm ${output_alignment}.unsorted
mv ${output_alignment}.sorted.bam ${output_alignment}

# mark duplicates
#java -jar /home/anduril/picard-tools-1.130/picard.jar MarkDuplicates \
#	INPUT=${output_alignment}.sorted \
#	OUTPUT=${output_alignment}.dupmarked \
#	METRICS_FILE=/home/anduril/test/test.picard.duplicate_metrics \
#	VALIDATION_STRINGENCY=LENIENT

# index
/home/anduril/samtools-1.2/samtools index ${output_alignment} 

# Archive log files
#mv *.out "${output_log}"

