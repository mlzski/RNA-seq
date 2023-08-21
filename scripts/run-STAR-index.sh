# script to run STAR to create a genome index for the aligner
# last updated: 2023-08-18

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: \
    (1) Path to data folder where index file will be placed, \
    (2) path to a .txt file that includes a list of genome reference files [.fa] and \
    (3) path to annotation file [.gtf] \
    ...Exiting"
    exit 1
fi	

# define arguments
index_dir=$1            # /path/to/index
ref_file=$2		# /path/to/reference_genome.fa
anno_file=$3            # /path/to/annotations.gtf

#NOTE: if multiple genome files (e.g., one per chromosome)
#ref_file=$(cat $2)      # /path/to/genome/fasta1 /path/to/genome/fasta2 ...

# export software (STAR)
export PATH=/nobackup/ummz/tools/bioinfo/STAR-2.7.11a/source:$PATH

# run STAR to generate index file
STAR \
   --runThreadN 8 \
   --runMode genomeGenerate \
   --genomeDir $index_dir \
   --genomeFastaFiles $ref_file \
   --sjdbGTFfile $anno_file \
   --sjdbOverhang 75

