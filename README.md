## ddRAD-seq
A repository for performing double-digest Restriction Associated DNA sequencing (ddRAD-seq) for reduced representation sequencing.

## Background
Double digest restriction-site associated DNA sequencing (ddRAD-seq) involves the use of two restriction enzymes to digest genomic DNA and create fragments of different sizes. The fragments are then size-selected to include only those in a specific length range. This selective process ensures that only a subset of the genome is sequences which reduces complexity and cost. Adapters are ligated to the size-selected fragments such that there is a different adapter for the overhang produced by the respective restriction enzyme. These are then amplified by PCR and sequenced using next-generation sequencing technology.

This method allows for reduced representation sequencing, which is a cost-effective alternative to whole genome sequencing. By using ddRAD-seq, we can identify single nucleotide polymorphisms (SNPs) without sequencing the entire genome. This technique provides control over the fraction of the genome sequenced, making it suitable for lower throughput applications. This approach allows for strain-level tracking and population structure analysis, achieving the necessary genetic insights within certain budget and technical constraints.

![image](https://github.com/KPU-AGC/ddRAD-seq/assets/90236200/953463f5-b1b8-4bf5-8833-319584008db7)

**Figure.** Diagram of ddRAD-seq workflow showing the construction of adapter ligation on double-digested fragments[^1].

## Objective
To simulate double digestion to identify the optimal pair of enzymes for performing ddRAD-seq on a given reference genome or simulated genome given genome length and GC content.

## Requirements
The recommended way to install dependencies is via `conda` using the `environment.yml` provided in the root directory of this repository.

The following programs and packages should be installed.
- `julia`
- `parallel`
- `fastx`
- `argparse`

## Usage
### 1A. Generating restriction enzyme statistics from a reference genome.
This method is probably the most reliable. With a very solid reference, this provides consistent cut sites and accurate prediction of restriction fragment loci.
1. (recommended) Put the reference genome in the directory provided.
2. Run the following command to generate statistics for a single run of double digestion of specific enzymes. 
```bash
REF="/path/to/reference_genome.fa"
ENZYME_1="EcoRI"
ENZYME_2="MseI"
julia src/jl-simRAD.jl use-ref ${REF} ${ENZYME_1} ${ENZYME_2}

# Or with size-selection (300 bp < x < 600 bp):
julia src/jl-simRAD.jl use-ref ${REF} ${ENZYME_1} ${ENZYME_2} -m 300 -M 600
```

### 1B. Generating restriction enzyme statistics without a reference genome.
If no reference genome is available, this program has the capability to simulate a genome with a given length and GC content. Compared to a stable reference genome, the results of this method are very inconsistent and should be run multiple times to be more accurate to the mean number of fragments produced.

1. Run the following command to generate statistics for a single run of double digestion of specific enzymes, similar to above.
```bash
GENOME_SIZE=46600000
GC_CONTENT=46
ENZYME_1="EcoRI"
ENZYME_2="MseI"
julia src/jl-simRAD.jl no-ref ${GENOME_SIZE} ${GC_CONTENT} ${ENZYME_1} ${ENZYME_2}

# Or with size-selection (300 bp < x < 600 bp):
julia src/jl-simRAD.jl no-ref ${GENOME_SIZE} ${GC_CONTENT} ${ENZYME_1} ${ENZYME_2} -m 300 -M 600
```

### 2. Running simulations of enzyme combinations in parallel.
To efficiently accomplish the goal of checking enzyme combinations, this program should be run in parallel using `GNU parallel`. There is a script provided in `src/` which takes an input list of enzymes and outputs a list of all unique combinations without replacement.

```bash
# 1. For use with a reference genome:
REF="/path/to/reference_genome.fa"
ENZYMES_LIST="/path/to/enzymes_list.txt"
./src/generate_combinations.sh | parallel --colsep '\t' --progress \
    'julia src/jl-simRAD.jl use-ref '"${REF}"' {1} {2} > output/{1}-{2}.csv'

# 2. For use with a simulated genome of 46.6 Mb, GC content of 46%, and with 10 repetitions:
GENOME_SIZE=46600000
GC_CONTENT=46
for i in $(seq 1 10); do
    mkdir -p gen_output/${i}
    ./src/generate_combinations.sh | parallel --colsep '\t' --progress \
        'julia src/jl-simRAD.jl no-ref '"${GENOME_SIZE}"' '"${GC_CONTENT}"' {1} {2} > gen_output/'"${i}"'/{1}-{2}.csv'
done
```
### 3. Aggregating results.
By default, results are output in `.csv` format for easy processing with R and other programs.
```bash
# 1. For reference-guided digestions:
find output -name "*.csv" -print | xargs cat > results.csv
# or honestly cat output/*.csv also works 

# 2. For reference-free digestions with multiple repetitions:
find gen_output -name "*.csv" -print | xargs cat  > results.csv

# And use this to add a header to the output.
sed -i '1s/^/enzymes,n_bases_covered,n_fragments_post_filter,n_fragments_pre_filter,genome_coverage\n/' results.csv
```

### 4. Output in GFF3 format.
When performing restriction digest with a given reference, it may be useful to see where restriction digestion occurs relative to genomic elements like genes or CpG islands. This makes the most sense, I think once the respective combination of restriction fragments has already been chosen.
```bash
REF="/path/to/reference_genome.fa"
ENZYME_1="EcoRI"
ENZYME_2="MseI"
julia src/jl-simRAD.jl use-ref ${REF} ${ENZYME_1} ${ENZYME_2} --gff3 true --csv false --pretty false > results.gff3
```
 

## References
[^1]: Liu, Michael & Worden, Paul & Monahan, Leigh & DeMaere, Matthew & Burke, Catherine & Djordjevic, Steven & Charles, Ian & Darling, Aaron. (2017). Evaluation of ddRADseq for reduced representation metagenome sequencing. PeerJ. 5. e3837. 10.7717/peerj.3837.
