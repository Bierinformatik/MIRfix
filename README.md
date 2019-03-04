# MIRfix
The MIRfix pipeline

## DESCRIPTION
MIRfix automatically curates miRNA datasets by improving alignments of
their precursors, the consistency of the annotation of mature miR and
miR* sequence, and the phylogenetic coverage. MIRfix produces
alignments that are comparable across families and sets the stage for
improved homology search as well as quantitative analyses.

## DEPENDENCIES
The easiest way to install MIRfix or its dependencies is using bioconda.
For this to work please first install [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html),
then [configure conda](https://conda.io/projects/conda/en/latest/user-guide/configuration/use-condarc.html) to use the [bioconda channel](https://bioconda.github.io/).

In short, first install conda, then run


```conda config --add channels defaults```


```conda config --add channels bioconda```


```conda config --add channels conda-forge```


Once conda is installed and configured, you can simply install MIRfix:

```conda create -n mirfix mirfix ``` will install MIRfix in the conda environment mirfix.

You can than use ```conda activate mirfix``` and run ```bash runMIRfix.sh``` from there.


Or you install all dependencies for MIRfix into a conda environment as follows:

In the *envs* directory you can find a file *MIRfix.env*, a simple

```conda env create -f MIRfix.env``` will install all dependencies in the conda environment mirfix.

You can than use ```conda activate mirfix``` and run ```bash runMIRfix.sh``` from a local clone of the git repository or the unpacked [source](https://github.com/Bierinformatik/MIRfix/archive/v1.0.tar.gz).

## RUNNING MIRFIX

Make sure to fix parameters and add files according to what you want to analyze.

For a test, go to the test directory and run ```bash testMIRfix.sh``` with files and parameters for testing or the ones already provided.

## PARAMETERS
### runMIRfix.sh
Following parameters are handled by the bash wrapper around MIRfix.py:

-location: Working directory, the path to where input, output, family sequences, family names, mapping files, genomes and mature sequences are located

-cores: Number of cores for parallel processing

-extension: Number of nucleotides for extension/trimming of precursor from genomic sequence

### MIRfix.py
For users who want to run the pipeline directly, following Parameters need to be set:

(cores) (outdir) (familydir) (list of families) (genomes) (mapfile) (matfile) (flank) [optional:(matrices)]

-(cores): Number of cores used for multiprocessing

-(outdir): output directory

-(familydir): Directory to fasta files of families to analyze

-(list of families): The list of families to analyze

-(genomes): The paths to the genomes where flanking regions to miRNAs can be extracted from

-(mapfile): mapping file between miRNA and miR-families

-(matfile): The mature sequences fasta file

-(flank): Number of precursor flanking nucleotides to take into account, we are using 10 in our study

-(matrices): The directory where the dialign2 matrices are located, please add the "/" at the end as in the example below *Please note that when using conda there is no need for this parameter*

## USAGE
```bash runMIRfix.sh LOCATION CORES EXTENSION```
To run an example, go to the directory "examples" and run the pipeline with ```bash testMIRfix.sh```.
You can also find example files for a run there.

## MORE
Please have a look at the detailed documentation in the docs folder.
