#!/usr/bin/env bash

location=${1:-${PWD}}				# Working directory DEFAULT is current working directory ${PWD}
cores=${2:-1}					# Numer of cores used for multithreading DEFAULT is 1
extension=${3:-10}				# Number of nucleotides for extension/trimming of precursor from genomic sequence

echo "Creating output directory ${location}/output"
mkdir -p ${location}/output

echo "Running MIRfix with ${cores} cores, ${extension}nt extension at ${location}"

cp mature.fa maturetest.fa
cp mapping.txt mappingtest.txt

python2 MIRfix.py ${cores} ${location}/output/ ${location}/Families/ ${location}/list.txt ${location}/genomes_list.txt ${location}/mappingtest.txt ${location}/maturetest.fa ${extension}

echo "Found "`diff mature.fa maturetest.fa`
echo "At "`diff mapping.txt mappingtest.txt

echo -e "\nCleanup\n"
rm -f maturetest.fa mappingtest.txt`
