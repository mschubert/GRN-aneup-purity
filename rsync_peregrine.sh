#!/bin/bash

rsync -auvr --include='*.pdf' \
    --exclude '.snakemake' \
    --exclude 'networks/**.RData' \
    --include '*.RData' \
    --include '*.png' \
    --include '*.tsv' \
    --include='*/' \
    --include='docs/*' \
    --exclude='*' \
    rug:/data/p282396/GRN-aneup-purity/* $(dirname $0)

#rsync -avur figure/* yoda:/nfs/research2/saezrodriguez/mike/speed2/figure/
