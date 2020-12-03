#!/bin/bash

# Paths
sequences="./GISAID/sequences"
reference="./GISAID/wuhan_1_reference.txt"
tools_folder=./covid/tools

# Test
cd $sequences
for file in *.txt; do $tools_folder/EMBOSS-6.6.0/emboss/stretcher \
	-auto -stdout \
    -asequence $reference \
    -bsequence $file \
    -datafile EDNAFULL -gapopen 16 -gapextend 4 \
    -aformat3 fasta -snucleotide1 -snucleotide2 \
    -outfile /Users/theophileouellette/Research/covid/GISAID/alignments/$file;
done