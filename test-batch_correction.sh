#!/bin/bash

# Set paths
scriptdir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

# Plain test
Rscript $scriptdir/batch_correction_docker_wrapper.R --loess "TRUE" dataMatrix "$scriptdir/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$scriptdir/test-data/input-batchcorrection-sampleMetadata.tsv" variableMetadata "$scriptdir/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$scriptdir/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$scriptdir/output-batchcorrection-variableMetadata.tsv" graph_output "/tmp/test.pdf" rdata_output "/tmp/test.Rdata"   || exit 1

diff $scriptdir/output-batchcorrection-dataMatrix.tsv $scriptdir/test-data/output-batchcorrection-dataMatrix.tsv || exit 2

# Test with different sample type column name and tags
Rscript $scriptdir/batch_correction_docker_wrapper.R --loess "TRUE" dataMatrix "$scriptdir/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$scriptdir/test-data/input-batchcorrection-sampleMetadata-customSampleType.tsv" variableMetadata "$scriptdir/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$scriptdir/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$scriptdir/output-batchcorrection-variableMetadata.tsv" graph_output "/tmp/test.pdf" rdata_output "/tmp/test.Rdata" sample_type_col_name MySampType sample_type_tags blank=blanc,pool=lot,sample=echant   || exit 1

diff $scriptdir/output-batchcorrection-dataMatrix.tsv $scriptdir/test-data/output-batchcorrection-dataMatrix.tsv || exit 2
