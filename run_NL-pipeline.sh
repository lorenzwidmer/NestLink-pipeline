#!/usr/bin/env bash

nextflow run main.nf \
    -profile standard \
    -params-file kdelr.json \
    -resume \
    -with-report \