# RNA-Seq pipeline

## Description: 
RNA-Seq pipeline incorporating six different programs for sequencing data processing. It can be used for both single-end or paired-end data. The processing starts from raw FASTQ data with quality control check, using FastQC. The raw data is then trimmed using Trimmomatic, which is followed by another quality control check using FastQC. The trimmed data is then aligned to the reference genome using STAR and indexed by SAMtools, with duplicates marking by Picard Tools. The final step includes read quantification using the featureCounts function from the Subread package to obtain read counts. More information about the analysis steps and software used can be found in wiki pages.

## Installation: 

### Requirements

- Access to a Linux-based HPC service. (Tested on a system based on the CentOS 7 distribution)
- Reference genome. It can be downloaded from https://genome.ucsc.edu/index.html
- All programs need to be downloaded and installed prior to running the pipeline. (See wiki for more detail)

### Download

Clone this repository to your own directory

```
git clone https://github.com/mihaux/RNA-seq.git
```

## Usage: 

**>>>The next section is usage, in which you instruct other people on how to use your project after they’ve installed it. This would also be a good place to include screenshots of your project in action.<<<**

## Contributing: 

**>>>Larger projects often have sections on contributing to their project, in which contribution instructions are outlined. Sometimes, this is a separate file. If you have specific contribution preferences, explain them so that other developers know how to best contribute to your work. To learn more about how to help others contribute, check out the guide for setting guidelines for repository contributors.<<<**

## Credits:

**>>>Include a section for credits in order to highlight and link to the authors of your project.

## License: 

**>>>Finally, include a section for the license of your project. For more information on choosing a license, check out GitHub’s licensing guide!<<<**
