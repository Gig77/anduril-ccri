#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

READS=$( getinput reads )
echo "Reads are ${READS}"

OPTIONS=$( getparameter options )
TMP=$(tempfile)

DOCKER=$( getparameter docker )
if [ -n "$DOCKER" ]; then
    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
fi

set -ex

cd $( dirname ${output_alignment} )

# align
$DOCKER gsnap \
	--format=sam \
	--quality-protocol=sanger \
	--print-snps \
	--input-buffer-size=5000 \
	--use-splicing=g1k_v37.splicesites \
	--use-snps=g1k_v37.snp138 \
	${OPTIONS} \
	\<\(java -jar /home/anduril/picard-tools-1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${READS} FASTQ=/dev/stdout\) \
		\| /home/anduril/samtools-1.2/samtools view -Shb - \
		\| /home/anduril/samtools-1.2/samtools sort -@ 2 -m 2000000000 - ${TMP}.sorted

# try to avoid error about missing current directory...
pushd .. && popd 

# mark duplicates
$DOCKER java -jar /home/anduril/picard-tools-1.130/picard.jar MarkDuplicates \
        INPUT=${TMP}.sorted.bam \
        OUTPUT=${TMP}.dupmarked \
        METRICS_FILE=$(dirname ${output_alignment})/picard.duplicate_metrics \
        VALIDATION_STRINGENCY=LENIENT

rm ${TMP}.sorted.bam
mv ${TMP}.dupmarked ${output_alignment}

# index
$DOCKER /home/anduril/samtools-1.2/samtools index ${output_alignment} 
