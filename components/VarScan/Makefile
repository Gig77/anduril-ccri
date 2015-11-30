.PHONY: test
test:
	docker pull biowaste:5000/anduril/gsnap
	cd /tmp && anduril test GSNAP --auto-bundles -b ${ANDURIL_BUNDLES}/ccri -b ${ANDURIL_BUNDLES}/sequencing
