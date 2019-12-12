#!/usr/bin/env bash

location=${1:-${PWD}}				# Working directory DEFAULT is current working directory ${PWD}
cores=${2:-1}					# Numer of cores used for multithreading DEFAULT is 1
extension=${3:-10}				# Number of nucleotides for extension/trimming of precursor from genomic sequence

echo "Creating output directory ${location}/output"
mkdir -p ${location}/output

echo "Running MIRfix with ${cores} cores, ${extension}nt extension at ${location}"

python3 MIRfix.py ${cores} ${location}/output/ ${location}/Families/ ${location}/listoffamilies.txt ${location}/fasta_to_search.txt ${location}/mapping_between_precursor_and_families.txt ${location}/mature_sequences.fa ${extension}

