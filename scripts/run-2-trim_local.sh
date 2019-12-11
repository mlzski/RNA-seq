
################################################################################# 
# this script runs trimming (on raw fastq files) for all files (pairs) at once
#################################################################################

if [ $# != 2 ] ; then
    echo -e "ERROR: 2 arguments are required: (1) Path to data folder and (2) path to output folder ... Exiting"
    exit 1
fi	

# define arguments
data_dir=$1			# path to the folder with fastq files
out_dir=$2			# path the the folder where output will be placed

# export software (Trimmomatic)
#export PATH=/nobackup/ummz/tools/Trimmomatic-0.39:$PATH

# check if output folders exist and create them if needed
if [ ! -d "$out_dir/processed_fastq" ]; then
    echo -ne "\nCreating directory..."
    mkdir $out_dir/processed_fastq
    echo -e $out_dir/processed_fastq
else
    read -p "Output directory already exists! Overwrite it? (y|n)? " choice
    case "$choice" in
    y|Y ) echo -e "\nOverwrote output directory."; rm -r $out_dir/processed_fastq; mkdir $out_dir/processed_fastq; echo -e $out_dir/processed_fastq;;   
    n|N ) echo -e "\nExiting..."; exit 1;;
    * )   echo "invalid"; exit 1;;
    esac
fi


# get cores of Read1 and Read2 names (with full directory paths)
r1_core=$(ls $data_dir/*_R1_001.fastq.gz | head -n 1)
r2_core=$(ls $data_dir/*_R2_001.fastq.gz | head -n 1)

# get rid of directory path and '_0001.fastq.gz'
r1_core_mod=$(basename $r1_core _001.fastq.gz)
r2_core_mod=$(basename $r2_core _001.fastq.gz)

echo $r1_core

echo $r1_core_mod

# run the trimmomatic command for each pair of fastq files
for filename in $(ls $data_dir/*_R1_001.fastq.gz)
do
    echo $filename
    echo 
    #echo java -jar /nobackup/ummz/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 -phred33 $r1_core $r2_core $out_dir/processed_fastq/${r1_core_mod}_paired.fq $out_dir/processed_fastq/${r1_core_mod}_unpaired.fq $out_dir/processed_fastq/${r2_core_mod}_paired.fq $out_dir/processed_fastq/${r2_core_mod}_unpaired.fq
done

