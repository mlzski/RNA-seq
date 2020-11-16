# RNA-Seq pipeline

## Description: 
RNA-Seq pipeline incorporating six different programs for sequencing data processing. It can be used for both single-end or paired-end data. The processing starts from raw FASTQ data with quality control check, using FastQC. The raw data is then trimmed using Trimmomatic, which is followed by another quality control check using FastQC. The trimmed data is then aligned to the reference genome using STAR and indexed by SAMtools, with duplicates marking by Picard Tools. The final step includes read quantification using the featureCounts function from the Subread package to obtain read counts. More information about the analysis steps and software used can be found in wiki pages.

## Installation: 

### Requirements

- Access to a Linux-based HPC service. (Tested on a system based on the CentOS 7 distribution)
- Reference genome. It can be downloaded from https://genome.ucsc.edu/index.html
- All programs need to be downloaded and installed prior to running the pipeline. (See 'Software installation' section on wiki for more detail)

### Download

Clone this repository to your own directory

```
git clone https://github.com/mihaux/RNA-seq.git
```

## Usage: 

See 'Usage example' section on wiki for a detailed explanation.

## Contributing: 

The pipeline has been created and tested on ARC4, part of the High Performance Computing facilities at the University of Leeds, UK.

## License: 

**>>>Finally, include a section for the license of your project. For more information on choosing a license, check out GitHub’s licensing guide!<<<**
