#!/usr/bin/env bash

nextflow run main.nf \
    -profile local \
    -entry prepare_data \
    -resume \
    -with-report \
    -with-timeline \
    -ansi-log false