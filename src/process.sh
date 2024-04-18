#!/bin/bash

ref="/mnt/c/Users/erick/Documents/GitHub/ddRAD-seq/ref/GCF_017499595.1_MGC_Penvy_1_genomic.fna"

./src/generate_combinations.sh | parallel --colsep '\t' --progress \
    'julia src/jl-simRAD.jl use-ref '"$ref"' {1} {2} > output/{1}-{2}.csv'

for i in 1 2 3 4 5; do
    mkdir -p gen_output/${i}
    ./src/generate_combinations.sh | parallel --colsep '\t' --progress \
        'julia src/jl-simRAD.jl no-ref 46600000 46 {1} {2} > gen_output/'"${i}"'/{1}-{2}.csv'
done