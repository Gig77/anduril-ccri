#!/bin/bash

# running functions.sh also sets logfile and errorfile.
source "$ANDURIL_HOME/bash/functions.sh"
export_command

#DOCKER=$( getparameter docker )
LABEL=$( getparameter label )

#if [ -n "$DOCKER" ]; then
#    DOCKER="docker run -u anduril --volumes-from anduril-master $DOCKER"
#fi

set -ex

# copy input file into temporary directory because GSEA has problems with special characters in path, e.g. hyphens (-)
# it also insists that file extension is .rnk ...
TMPDIR=$(mktemp -d)
RNKFILE=$TMPDIR/$LABEL.rnk
cp ${input_rnk} $RNKFILE

cd $TMPDIR && java -cp /data_synology/software/gsea-2.0.13/gsea2-2.0.13.jar -Xmx3048m xtools.gsea.GseaPreranked \
		-rpt_label $LABEL \
		-rnk $RNKFILE \
		-gmx /mnt/projects/generic/data/msigdb5.0/msigdb.v5.0.symbols.gmt,/mnt/projects/iamp/results/gsea/custom/mir_target_genesets.gmx \
		-out $TMPDIR \
		-nperm 1000 \
		-plot_top_x 3000 \
		-collapse false \
		-mode Max_probe \
		-norm meandiv \
		-scoring_scheme weighted \
		-include_only_symbols true \
		-make_sets true \
		-rnd_seed 149 \
		-set_max 500 \
		-set_min 10 \
		-zip_report false \
		-gui false 

rename "s/$LABEL\\.GseaPreranked\\.\\d+\$/gsea_output/" $TMPDIR/$LABEL.GseaPreranked*
cp  $TMPDIR/gsea_output/gsea_report_for_na_pos_*.xls ${output_enrichedUp}.part
cp  $TMPDIR/gsea_output/gsea_report_for_na_neg_*.xls ${output_enrichedDown}.part
rm -rf $(dirname ${output_enrichedUp})/gsea_output
mv $TMPDIR/gsea_output $(dirname ${output_enrichedUp})/

Rscript ${metadata_componentPath}/add_genes.r \
	--gsea-result-file ${output_enrichedUp}.part \
	--phenotype up \
	--gene-set-dir $(dirname ${output_enrichedUp})/gsea_output \
	--num-genes $(cat ${input_rnk} | wc -l) \
	--output-file ${output_enrichedUp}
rm ${output_enrichedUp}.part

Rscript ${metadata_componentPath}/add_genes.r \
	--gsea-result-file ${output_enrichedDown}.part \
	--phenotype down \
	--gene-set-dir $(dirname ${output_enrichedUp})/gsea_output \
	--num-genes $(cat ${input_rnk} | wc -l) \
	--output-file ${output_enrichedDown}
rm ${output_enrichedDown}.part

rm -rf $TMPDIR