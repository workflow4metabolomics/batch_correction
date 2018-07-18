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

# Run MTBLS404 test {{{1
################################################################

run_mtbls404_test() {

	dir="$PROG_PATH/res/mtbls404"
	matrix_file="$dir/MTBLS404_W4M_data.tsv"
	samp_file="$dir/MTBLS404_W4M_samp.tsv"
	var_file="$dir/MTBLS404_W4M_var.tsv"
	matrix_outfile="$PROG_PATH/output-batchcorrection-dataMatrix.tsv"
	var_outfile="$PROG_PATH/output-batchcorrection-variableMetadata.tsv"
	graph_outfile="$PROG_PATH/output-graph.pdf"
	rdata_outfile="$PROG_PATH/output-rdata.Rdata"

	Rscript $ROOTDIR/batch_correction_all_loess_wrapper.R dataMatrix "$matrix_file" sampleMetadata "$samp_file" variableMetadata "$var_file" method "all_loess_pool" span "1.0" dataMatrix_out "$matrix_outfile" variableMetadata_out "$var_outfile" graph_output "$graph_outfile" rdata_output "$rdata_outfile" batch_col_name "Factor.Value.Batch." injection_order_col_name "Factor.Value.Injection.order." sample_type_col_name "Factor.Value.Material.type." || exit 1
}

# MAIN {{{1
################################################################

run_batch_correction_simple_test
run_same_test_as_batch_correction_XML
run_mtbls404_test
