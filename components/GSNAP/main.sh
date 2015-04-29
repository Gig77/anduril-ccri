#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"

export_command

set -ex

mkdir "${output_log}"

cd $( dirname ${output_alignment} )

READS=$( getinput reads )
echo "Reads are ${READS}"

OPTIONS=$( getparameter options )
CMD_PREFIX=$( getparameter cmdPrefix )

echo cmdPrefix: ${CMD_PREFIX}
 
# align
${CMD_PREFIX} "gsnap \
	--db=g1k_v37_etv6runx1 \
	--dir=/home/anduril/g1k_v37_etv6runx1 \
	--format=sam \
	--quality-protocol=sanger \
	--print-snps \
	--input-buffer-size=5000 \
	--use-splicing=g1k_v37.splicesites \
	--use-snps=g1k_v37.snp138 \
	${OPTIONS} \
	<(java -jar /home/anduril/picard-tools-1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${READS} FASTQ=/dev/stdout) \
		| /home/anduril/samtools-1.2/samtools view -Shb - \
		| /home/anduril/samtools-1.2/samtools sort -@ 2 -m 2000000000 - ${output_alignment}.sorted"

# mark duplicates
${CMD_PREFIX} java -jar /home/anduril/picard-tools-1.130/picard.jar MarkDuplicates \
        INPUT=${output_alignment}.sorted.bam \
        OUTPUT=${output_alignment}.dupmarked \
        METRICS_FILE=picard.duplicate_metrics \
        VALIDATION_STRINGENCY=LENIENT

rm ${output_alignment}.sorted.bam
mv ${output_alignment}.dupmarked ${output_alignment}

# index
${CMD_PREFIX} /home/anduril/samtools-1.2/samtools index ${output_alignment} 

# Archive log files
#mv *.out "${output_log}"

