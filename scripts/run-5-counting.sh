###################################################################################################### 
# this script runs Cufflinks to count reads from STAR output (name_Aligned.sortedByCoord.out.bam)
###################################################################################################### 

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: (1) path to data folder, (2) path to output directory, (3) path to index file and (4) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi    

# export software (Cufflinks)
export PATH=/nobackup/ummz/tools/cufflinks-2.2.1.Linux_x86_64:$PATH

# assign variables
data_dir=$1
out_dir=$2
index_dir=$3

# get the read1 fastq.gz file and its pair
fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")

read1=$fastqFile
read2=$(echo $read1 | sed 's/R1/R2/g')

bam_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

# run Cufflinks

cufflinks -p 4 -G /nobackup/<yourFolder>/data/annotation/human_chr22.gtf -o /nobackup data/cufflinks/normal_rep1_clout_with_G
/nobackup     data/bam/normal_rep1_Aligned.sortedByCoord.out.bam


