#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

QUERY=$( getinput query )
REFERENCE=$( getinput reference )
PVALUE=$( getparameter pvalue )
MAX_REGIONS=$( getparameter maxRegionsPerChr )
TMPDIR=$( gettempdir )

DOCKER=$( getparameter docker )
if [ -n "$DOCKER" ]; then
    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
fi

# NOTE: $(eval echo $CMD_PREFIX) substitutes the working directory variable in the prefix with the actual value; otherwise docker will complain

set -ex

# prepare latex output

rm -f ${output_document}/document.tex
mkdir -p ${output_document}
SECTION_TITLE=$( getparameter sectionTitle )
SECTION_TYPE=$( getparameter sectionType )
CAPTION=$( getparameter caption )

if [ -n "$SECTION_TYPE" ]; then
	SECTION_TYPE="\\$SECTION_TYPE"
fi

if [ -n "$SECTION_TITLE" ]; then
	echo "$SECTION_TYPE{$SECTION_TITLE}" | sed 's/\\\\/\\/g' >> ${output_document}/document.tex
fi

# perform enrichment (if input file is not empty)

cat ${QUERY} | sed '1d' | sed 's/\"//g' > ${TMPDIR}/query.txt

if [ -s "${TMPDIR}/query.txt" ]; then
	zcat -f ${REFERENCE} | perl -ne 'print "$4\t$1\t$2\t$3\n" if (/^([^\t]+).*\tgene\t(\d+)\t(\d+)\t.*gene_id "([^"]+)/)' > ${TMPDIR}/reference.txt
	
	$DOCKER perl /home/anduril/pge.pl -g -a ${PVALUE} -r user -f ${TMPDIR}/reference.txt -q ${TMPDIR}/query.txt > ${output_enrichedRegions}.part
	
	cat ${output_enrichedRegions}.part | (read -r; printf "%s\n" "$REPLY"; sort -k 4g) > ${output_enrichedRegions} # sort by p-value ascending
	rm ${output_enrichedRegions}.part
	
	# ideogram plot

	PDF=${output_document}/ideogram_${metadata_instanceName}.pdf
	$DOCKER Rscript /home/anduril/ideogram.r --query ${QUERY} --reference ${TMPDIR}/reference.txt --regions ${output_enrichedRegions} --output ${PDF} --maxRegionsPerChr ${MAX_REGIONS}

	# add latex fragment for figure
			
	echo '\begin{figure}[!ht]' >> ${output_document}/document.tex
	echo '\begin{center}' >> ${output_document}/document.tex
	echo '\includegraphics{'$(basename $PDF)'}' >> ${output_document}/document.tex
	echo '\end{center}' >> ${output_document}/document.tex
	echo "\caption{${CAPTION}}" | sed 's/\\\\/\\/g' >> ${output_document}/document.tex
	echo '\end{figure}' >> ${output_document}/document.tex
else
	echo 'No DEGs found.' >>  ${output_document}/document.tex
	echo -e "chr\tstart\tend\tpvalue\tpvalueadj\tcommon\tsize\tgenes" > ${output_enrichedRegions}
fi


