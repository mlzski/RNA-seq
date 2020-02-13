###################################################################################################### 
# this script runs Cufflinks to count reads from STAR output (name_Aligned.sortedByCoord.out.bam)
###################################################################################################### 

if [ $# != 4 ] ; then
    echo -e "ERROR: 4 arguments are required: (1) path to data folder, (2) path to output directory, (3) path to annotation (.gtf) file and (4) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi    

# export software (Cufflinks)
export PATH=/home/home02/ummz/tools/cufflinks-2.2.1.Linux_x86_64:$PATH

# assign variables
data_dir=$1
out_dir=$2
anno_file=$3

# get the input (bam file) names
bam_file=$(ls $data_dir/*Aligned.sortedByCoord.out.bam | sed -n -e "$SGE_TASK_ID p")

# get the 'core' of file names
core_name=$(ls $data_dir/*_Aligned.sortedByCoord.out.bam | rev | cut -d '/' -f 1 | cut -c 31- | rev | sed -n -e "$SGE_TASK_ID p")

# create directories for output files (one pre sample)
mkdir $out_dir/out_${core_name}

# run Cufflinks
cufflinks -p 8 -G $anno_file -o $out_dir/out_${core_name} $bam_file

