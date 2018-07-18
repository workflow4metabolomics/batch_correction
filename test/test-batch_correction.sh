#!/bin/bash

# Constants {{{1
################################################################

PROG_PATH=$(dirname $0)
ROOTDIR="$PROG_PATH/.."

# Run same test as Batch correction XML {{{1
################################################################

run_same_test_as_batch_correction_XML() {

	Rscript $ROOTDIR/batch_correction_all_loess_wrapper.R dataMatrix "$ROOTDIR/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$ROOTDIR/test-data/input-batchcorrection-sampleMetadata.tsv" variableMetadata "$ROOTDIR/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$PROG_PATH/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$PROG_PATH/output-batchcorrection-variableMetadata.tsv" graph_output "$PROG_PATH/output-graph.pdf" rdata_output "$PROG_PATH/output-rdata.Rdata" batch_col_name batch injection_order_col_name injectionOrder sample_type_col_name sampleType || exit 1

	diff "$PROG_PATH/output-batchcorrection-dataMatrix.tsv" "$ROOTDIR/test-data/output-batchcorrection-dataMatrix.tsv" || exit 2
}

# Run batch correction simple test {{{1
################################################################

run_batch_correction_simple_test() {
	Rscript $ROOTDIR/batch_correction_docker_wrapper.R --loess "TRUE" dataMatrix "$ROOTDIR/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$ROOTDIR/test-data/input-batchcorrection-sampleMetadata.tsv" variableMetadata "$ROOTDIR/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$PROG_PATH/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$PROG_PATH/output-batchcorrection-variableMetadata.tsv" graph_output "$PROG_PATH/output-graph.pdf" rdata_output "$PROG_PATH/output-rdata.Rdata" || exit 1

	diff $PROG_PATH/output-batchcorrection-dataMatrix.tsv $ROOTDIR/test-data/output-batchcorrection-dataMatrix.tsv || exit 2
}

# MAIN {{{1
################################################################

run_batch_correction_simple_test
run_same_test_as_batch_correction_XML
