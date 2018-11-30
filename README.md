# MIRfix
The MIRfix pipeline

##DESCRIPTION
MIRfix automatically curates miRNA datasets by improving alignments of
their precursors, the consistency of the annotation of mature miR and
miR* sequence, and the phylogenetic coverage. MIRfix produces
alignments that are comparable across families and sets the stage for
improved homology search as well as quantitative analyses.

##DEPENDENCIES
The easiest way to install dependencies is using bioconda.
In the envs directory you can find a file MIRfix.env, a simple ```conda env create -f MIRfix.env``` will install all dependencies in the conda environment mirfix.
You can than use ```conda activate mirfix``` and run ```bash runMIRfix.sh``` from there.
Make sure to fix parameters and add files according to what you want to analyze.
For a test, go to the test directory and run ```bash testMIRfix.sh``` with files and parameters for testing or the ones already provided.

##PARAMETERS
(1 or 2) (outdir) (filesdir) (file or list) (genomes) (mapfile) (matfile) (matrice) (flank)

(1 or 2): 1 for a file and for list of files
(outdir): output directory
(filesdir): In case of list, the directory where all the files are, in case of file just write "null"
(file or list): Just add the list in case of list and file in case of only file
(genomes): The genomes to search (list) as explained in the Documentation
(mapfile): mapping file
(matfile): The mature sequences fasta file
(matrices): The directory where the matrices are, please add the "/" at the end as in the example below
(flank): we are using 10 in our study

##USAGE
bash runMIRfix.sh
