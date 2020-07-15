
######################################################################################################
# this script runs STAR to create a genome index for the aligner
######################################################################################################

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: \
    (1) Path to data folder where index file will be placed, \
    (2) path to genome reference file [.fa] and \
    (3) path to annotation file [.gtf] \
    ...Exiting"
    exit 1
fi	

# define arguments
index_dir=$1            # /path/to/genomeDir
ref_file=$2             # /path/to/genome/fasta1 /path/to/genome/fasta2 ...
anno_file=$3            # /path/to/annotations.gtf

# export software (STAR)
export PATH=/home/home02/ummz/tools/STAR-2.7.3a/bin/Linux_x86_64_static:$PATH

# run STAR to generate index file
STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir $index_dir \
--genomeFastaFiles $ref_file \
--sjdbGTFfile $anno_file \
--sjdbOverhang 150

# descriptions:
# --runThreadN          => defines the number of threads to be used for genome generation, it has to be set to the number of available cores on the server node.
# --runMode             => genomeGenerate option directs STAR to run genome indices generation job.
# --genomeDir           => specifies path to the directory (henceforth called ”genome directory” where the genome indices are stored. This directory has to be created (with mkdir) before STAR run and needs to have writing permissions. The file system needs to have at least 100GB of disk space available for a typical mammalian genome. It is recommended to remove all files from the genome directory before running the genome generation step. This directory path will have to be supplied at the mapping step to identify the reference genome.
# --genomeFastaFiles    => specifies one or more FASTA files with the genome reference sequences. Multiple reference sequences (henceforth called chromosomes) are allowed for each fasta file. You can rename the chromosomes names in the chrName.txt keeping the order of the chromo- somes in the file: the names from this file will be used in all output alignment files (such as .sam). The tabs are not allowed in chromosomes names, and spaces are not recommended.
# --sjdbGTFfile         => specifies the path to the file with annotated transcripts in the standard GTF format. STAR will extract splice junctions from this file and use them to greatly improve accuracy of the mapping. While this is optional, and STAR can be run without annotations, using annotations is highly recommended whenever they are available. Starting from 2.4.1a, the annotations can also be included on the fly at the mapping step.
# --sjdbOverhang        => specifies the length of the genomic sequence around the annotated junction to be used in constructing the splice junctions database. Ideally, this length should be equal to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina 2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as well as the ideal value.
