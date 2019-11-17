#!/usr/bin/env bash

#test="RF00725_MI0001598"
test="RF00725_MI0001636"
location=${1:-${PWD}}				# Working directory DEFAULT is current working directory ${PWD}
cores=${2:-1}					# Numer of cores used for multithreading DEFAULT is 1
extension=${3:-10}				# Number of nucleotides for extension/trimming of precursor from genomic sequence

echo "Running MIRfix with ${cores} cores, ${extension}nt extension at ${location}"

cp $location/${test}_mature_original.fa $location/${test}_maturetest.fa
cp $location/${test}.mapping_original.txt $location/${test}.mappingtest.txt

python3 MIRfix.py3 -j ${cores} -o ${location}/Test_Cristian_output/ -i ${location}/Families/ -f ${location}/${test}_list.txt -g ${location}/${test}_genomes_list.txt -m ${location}/${test}.mappingtest.txt -a ${location}/${test}_maturetest.fa -e ${extension}

#python2 MIRfix.py ${cores} ${location}/output/ ${location}/Families/ ${location}/${test}_list.txt ${location}/${test}_genomes_list.txt ${location}/${test}.mappingtest.txt ${location}/${test}_maturetest.fa ${extension}

#echo "Found "`diff mature.fa maturetest.fa`
#echo "At "`diff mapping.txt mappingtest.txt

#echo -e "\nCleanup\n"
#rm -f maturetest.fa mappingtest.txt`
