### this script runs Salmon to obtaine reads from fastq files (trimmed) on transcript level
# NOTICE: prior to launching this command, make sure that the index has been already generated (just with a single command) 

# TODO: finish implementation for single-end data

# this version is based on this tutorial: https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/salmon_kalliso_STAR_compare.md

if [ $# != 6 ] ; then
    echo -e "ERROR: 6 arguments are required: \
    (1) running mode (either 'SE' or 'PE'), \
    (2) path to folder with input data, \
    (3) path to folder for output directory, \
    (4) path to index file \
    (5) path to .gtf file and \
    (6) SGE_TASK_ID argument for array jobs ...Exiting"
    exit 1
fi

# export software (salmon)
#export PATH=/home/home02/ummz/tools/bioinfo/salmon-latest_linux_x86_64/bin:$PATH
export PATH=/nobackup/ummz/tools/bioinfo/salmon-1.6.0_linux_x86_64/bin:$PATH

# assign variables
run_mode=$1               # 'transcript-level' or 'gene-level'
data_dir=$2
out_dir=$3
index_dir=$4
gtf_dir=$5

# prepare params.txt filr with information about this run
me=$(basename "$0")
echo -ne "Script directory:" `pwd`"/"$me "\n" >> params_${me}_.txt

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`";
echo -ne "Executed on:" $current_date_time "\n" >> params_${me}_.txt

# run Salmon mapping job
if [ $run_mode == 'transcript-level' ] ; then
   
   # get the read1 fastq.gz file and its pair read2
    fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
    samp_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

    # run Salmon in paired-end mode [PE]
    salmon quant \
    	-i $index_dir \
	-l ISR \
	-1 $read1 \
	-2 $read2 \
	-o $out_dir/ID-${samp_name} \
	-p 8
    
    # NOTE: use either "-l ISR" or "-l A" (automatic); for new data, always "A"
 
    # print the main command to params.txt
    echo "salmon quant -i $index_dir -l ISR -1 $read1 -2 $read2 -p 8 -o $out_dir/${samp_name}_quant" >> params_${me}_.txt


elif [ $run_mode == 'gene-level' ] ; then

    # get the read1 fastq.gz file and its pair read2
    fastqFile=$(ls $data_dir/*_R1_paired.fq | sed -n -e "$SGE_TASK_ID p")
    read1=$fastqFile
    read2=$(echo $read1 | sed 's/R1/R2/g')
    samp_name=$(ls $data_dir/*_R1_paired.fq | rev | cut -d '/' -f 1 | cut -c 13- | rev | sed -n -e "$SGE_TASK_ID p")

    # run Salmon in paired-end mode [PE]
    salmon quant \
	-i $index_dir \
	-l ISR \
	-1 $read1 \
	-2 $read2 \
	-o $out_dir/ID-${samp_name} \
	-g $gtf_dir \
	--seqBias \
	--validateMappings \
	-p 8

    # NOTE: use either "-l ISR" or "-l A" (automatic); for new data, always "A"
 
    # print the main command to params.txt
    echo "salmon quant -i $index_dir -l ISR -1 $read1 -2 $read2 -o $out_dir/ID-${samp_name} -g $gtf_dir --seqBias --validateMappings -p 8" >> params_${me}_.txt
      

# NOTES:
# --seqBias 
# --validateMappings

else
    echo "ERROR... run_mode argument must be specified as: 'transcript-level' or 'gene-level'"
    exit 1
fi

echo "Finished for " $samp_name 
