#!/bin/bash

# Set paths
scriptdir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

Rscript $scriptdir/batch_correction_docker_wrapper.R --loess "TRUE" dataMatrix "$scriptdir/test-data/input-batchcorrection-dataMatrix.tsv" sampleMetadata "$scriptdir/test-data/input-batchcorrection-sampleMetadata.tsv" variableMetadata "$scriptdir/test-data/input-batchcorrection-variableMetadata.tsv" method "all_loess_pool" span "1" dataMatrix_out "$scriptdir/output-batchcorrection-dataMatrix.tsv" variableMetadata_out "$scriptdir/output-batchcorrection-variableMetadata.tsv" graph_output "/tmp/test.pdf" rdata_output "/tmp/test.Rdata"   || exit 1

diff $scriptdir/output-batchcorrection-dataMatrix.tsv $scriptdir/test-data/output-batchcorrection-dataMatrix.tsv || exit 2
