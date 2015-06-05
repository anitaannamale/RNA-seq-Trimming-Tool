# RNA-seq-Trimming-Tool

## About
Most modern sequencing technologies produce read that have a deteriorating quality towards the ends, incorrectly called bases in some regions and ... . All these bad quality reads negatively impact on assembly, mapping and downstream bioinfomatics analyses. These are particulary true, in de novo assembly where the genome or the transcriptome is reconstructed only based on information contained in reads.

This is an independant module for adapter and quality trimming of RNA-seq data which uses Trimmomatic. This module is integrated in a pipeline of de novo assembly for non-models organisms. Here is the pipeline link : https://github.com/arnaudmeng/denovo-assembly-pipeline-upmc

## Version
0.0.1

## Requirements

Python 2.7.6 & Java (Trimmomatic)

## Usage

This module has two ways of working (reading input files and trimming parameters) from : 
      - an XML file
      - or command line directly

And each way has two modes to work :
     - single ends "SE" `./Filtrage SE`
     - paired ends "PE" `./Filtrage PE`

To read input files and trimming parameters from the XML file, run :

            'python ./Filtrage.py --XML'

To launch trimmomatic using commandline arguments see below:

For more informations, see module help, running:

            'python ./Filtrage.py -h'
            
## Single ends data

`./Filtrage.py SE` takes an input fastq file and outputs a trimmed version of that file. It has all options given by Trimmomatic :



## Paired ends data

