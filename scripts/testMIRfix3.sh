#!/usr/bin/env bash

location=${1:-${PWD}}				# Working directory DEFAULT is current working directory ${PWD}
cores=${2:-1}					# Numer of cores used for multithreading DEFAULT is 1
extension=${3:-10}				# Number of nucleotides for extension/trimming of precursor from genomic sequence

echo "Running MIRfix with ${cores} cores, ${extension}nt extension at ${location}"

cp mature.fa maturetest.fa
cp mapping.txt mappingtest.txt

python3 MIRfix.py3 -j ${cores} -o ${location}/Test_output/ -i ${location}/Families/ -f ${location}/list.txt -g ${location}/genomes_list.txt -m ${location}/mappingtest.txt -a ${location}/maturetest.fa -e ${extension}

echo "Found "`diff mature.fa maturetest.fa`
echo "At "`diff mapping.txt mappingtest.txt

echo -e "\nCleanup\n"
rm -f maturetest.fa mappingtest.txt`
