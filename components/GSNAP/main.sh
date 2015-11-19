#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

READS=$( getinput reads )
echo "Reads are ${READS}"

EXECUTABLE=$( getparameter executable )
OPTIONS=$( getparameter options )
INPUTTYPE=$( getparameter inputType )
TMP=$(tempfile)

DOCKER=$( getparameter docker )
if [ -n "$DOCKER" ]; then
    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
fi

set -ex

cd $( dirname ${output_alignment} )


# deduce input file type from file name suffix
if [ "${INPUTTYPE}" = "auto" ]; then 
	if [ "${READS#*.}" = "bam" ]; then
		INPUTTYPE='BAM'
	elif [ "${READS#*.}" = "fastq" ] || [ "${READS#*.}" = "fastq.gz" || [ "${READS#*.}" = "txt.gz" ]; then
		INPUTTYPE='FASTQ'
	else
		echo "ERROR: Could not deduce input file type from file name: ${READS}"
		exit 1
	fi
fi

# align
if [ "${INPUTTYPE}" = "BAM" ]; then 
	$DOCKER $EXECUTABLE \
		--format=sam \
		--quality-protocol=sanger \
		--print-snps \
		--input-buffer-size=5000 \
		${OPTIONS} \
		\<\(java -jar /home/anduril/picard-tools-1.130/picard.jar SamToFastq VALIDATION_STRINGENCY=SILENT INPUT=${READS} FASTQ=/dev/stdout\) \
			\| /home/anduril/samtools-1.2/samtools view -Shb - \
			\| /home/anduril/samtools-1.2/samtools sort -@ 2 -m 2000000000 - ${TMP}.sorted
else
	$DOCKER $EXECUTABLE \
		--format=sam \
		--quality-protocol=sanger \
		--print-snps \
		--input-buffer-size=5000 \
		${OPTIONS} \
		\<\(gunzip -fc ${READS}\) \
			\| /home/anduril/samtools-1.2/samtools view -Shb - \
			\| /home/anduril/samtools-1.2/samtools sort -@ 2 -m 2000000000 - ${TMP}.sorted
fi

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
