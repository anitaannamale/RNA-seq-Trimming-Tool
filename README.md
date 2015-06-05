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

- remove adapter sequences : `illuminaclip <fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>`
- quality trimming : `-slidingwindow <window-size>:<required-quality>`
- adaptative quality trimming depending on the length of reads : `-maxinfo <target-length>:<strictness>`
- trim base from 5' end until the required minimal quality is achieved : `-leading <required-quality>`
- trim base from 3' end until the required minimal quality is achieved : `-trailing <required-quality>`
- trim a fixed number of bases from 5' end : `-head <number>`
- trim a fixed number of bases from 3' end : `-crop <number>`
- remove read shorter than a given length : `-minlen <length>`

### Examples :

      python Filtrage.py SE read_1.fastq -illuminaclip fasta-file.fa:2:10:30
      python Filtrage.py SE read_1.fq -slidingwindow 10:30 -leading 30 -minlen 36
      
## Paired ends data

`./Filtrage.py PE` takes two input fastq files and outputs a trimmed version of these files and a two files containing 'single reads' for each direction. 'Single reads' are reads who passed the filter for one direction but not the other. This mode has also all options given by Trimmomatic, see above.

### Examples :

      python Filtrage.py PE read_1.fastq read_2.fastq -illuminaclip fasta-file.fa:2:10:30 -crop 10
      python Filtrage.py PE read_1.fq.bz2 read_2.fq.bz2 -illuminaclip fasta-file.fa:2:10:30 -slidingwindow 10:30 -minlen 36
