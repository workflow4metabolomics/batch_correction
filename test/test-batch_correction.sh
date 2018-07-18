#!/bin/bash

# Constants {{{1
################################################################

PROG_PATH=$(dirname $0)
ROOTDIR="$PROG_PATH/.."

# Run same test as Galaxy XML {{{1
################################################################

run_same_test_as_Galaxy_XML() {
	Rscript $ROOTDIR/batch_correction_docker_wrapper.R --loess "TRUE" dataMatrix "$ROOTDIR/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$ROOTDIR/test-data/input-batchcorrection-sampleMetadata.tsv" variableMetadata "$ROOTDIR/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$PROG_PATH/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$PROG_PATH/output-batchcorrection-variableMetadata.tsv" graph_output "/tmp/test.pdf" rdata_output "/tmp/test.Rdata"   || exit 1

	diff $PROG_PATH/output-batchcorrection-dataMatrix.tsv $ROOTDIR/test-data/output-batchcorrection-dataMatrix.tsv || exit 2
}

# MAIN {{{1
################################################################

run_same_test_as_Galaxy_XML
