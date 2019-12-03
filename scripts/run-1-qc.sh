
################################################################################# 
# this script runs first QC analysis (on raw fatsq files) for all files at once
#################################################################################

if [ $# != 2 ] ; then
    echo -e "ERROR: 2 arguments are required: (1) Path to data folder and (2) path to output folder ... Exiting"
    exit 1
fi	

# TODO: add count to see how many fastq files are in the directoey and to track their processing


# define arguments
data_dir=$1			# path to the folder with fastq files
out_dir=$2			# path the the folder where output will be placed

# export software (FastQC)
export PATH=/nobackup/leedsomics_workshop/tools/FastQC:$PATH

# check if output folders exist and create them if needed
if [ ! -d "$out_dir/report" ]; then
    echo -ne "\nCreating directories..."
    mkdir $out_dir/report $out_dir/temp
    echo -e $out_dir/report
    echo -e $out_dir/temp
else
    read -p "Output directory already exists! Overwrite it? (y|n)? " choice
    case "$choice" in
    y|Y ) echo -e "\nOverwrote output directory."; rm -r $out_dir/report $out_dir/temp; mkdir $out_dir/report $out_dir/temp; echo -e $out_dir/report; echo -e $out_dir/temp;;   
    n|N ) echo -e "\nExiting..."; exit 1;;
    * )   echo "invalid"; exit 1;;
    esac
fi
 

for filename in $(ls $data_dir/*.fastq.gz) 
do
    fastqc -o $out_dir/report --threads 4 --dir $out_dir/temp $filename
done











