
######################################################################################################
# this script runs STAR to create a genome index for the aligner
######################################################################################################

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: (1) Path to data folder where index file will be placed, (2) path to genome reference file [.fa] and (3) path to annotation file [.gtf] ... Exiting"
    exit 1
fi	

# define arguments
index_dir=$1
ref_file=$2
anno_file=$3

# export software (STAR)
export PATH=/nobackup/ummz/tools/STAR-2.7.3a/bin/Linux_x86_64_static:$PATH

# run STAR to generate index file
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index_dir --genomeFastaFiles $ref_file --sjdbGTFfile $anno_file --sjdbOverhang 150
