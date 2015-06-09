.PHONY: test
test:
	docker pull biowaste:5000/anduril/rseqc
	cd /tmp && anduril test CalcGeneBodyCoverage --auto-bundles -b ${ANDURIL_BUNDLES}/ccri -b ${ANDURIL_BUNDLES}/sequencing
