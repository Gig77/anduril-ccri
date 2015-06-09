.PHONY: test
test:
	docker pull biowaste:5000/anduril/pge
	cd /tmp && anduril test PGE --auto-bundles -b ${ANDURIL_BUNDLES}/ccri -b ${ANDURIL_BUNDLES}/sequencing
