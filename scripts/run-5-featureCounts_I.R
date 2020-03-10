###################################################################################################### 
# this script runs featureCounts to count reads from STAR output (name_Aligned.sortedByCoord.out.bam)
# using the array concept with the 'SGE_TASK_ID' parameter
# so it launches all files (samples) in parallel and produces one output (.csv) file per each sample 
###################################################################################################### 

# INPUT:
# (1) input a set of files that contain read mapping results
# (2) an annotation file that includes genomic features. (either in SAM and BAM format, which is detected automatically) 
# NOTICE: Input reads can be name-sorted or location-sorted. 

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!requireNamespace("Rsubread", quietly = TRUE)) BiocManager::install("Rsubread")
library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)

# for running in R-studio
# args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam", 
#          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/5_/postprocessed")

if (length(args)!=4) {
  stop("4 arguments must be supplied: 
	\n(1) running mode parameter [SE or PE],
	\n(2 - input) path to directory with data, 
	\n(3 - output) path where output files should be stored and 
	\n(4) SGE_TASK_ID argument for array jobs", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[2], sep="\n")

cat("Directory for results (OUT): ")
cat(args[3], sep="\n")

setwd(args[3])

# get task ID
task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# retrive built-in annotation file (other genomes available: mm9, mm10, hg19 and hg38 (NCBI RefSeq))
ann <- getInBuiltAnnotation(annotation = "hg38")

# create a list with file names to be processed:
files_list <- list.files(args[2], pattern = ".bam$", full.names = TRUE)

# create the parameter for running as an array
parameter <- files_list[task_id]

# run featureCounts() command [SE - single end | PE - paired end]

if(args[1] == "SE"){
	fc_SE <- featureCounts(parameter, annot.ext=ann, nthreads = 40)	
} else if (args[1] == "PE"){
	fc_PE <- featureCounts(parameter, annot.ext=ann, isPairedEnd=TRUE, nthreads = 40)	
} else {
	stop("ERROR: running mode parameter must be defined as either SE or PE", call.=FALSE)	
}

# write output table with counts and statistics
outputName_1 <- paste("counts_",task_id,".csv",sep="")
write.csv(fc_SE$counts, file=outputName_1)

outputName_2 <- paste("stats_",task_id,".csv",sep="")
write.csv(fc_SE$stat, file=outputName_2)


######################################################################################################
# featureCOunts() OUTPUTS:
# fc_SE$counts      => a very long vector with counts for all genes
# fc_SE$annotation  => a table with all genes (rows) 6 columns for: "GeneID", "Chr", "Start", "End", "Strand", "Length"
# fc_SE$targets     => file name(s) used
# fc_SE$stat        => a table with statistics (14 metrics)

# return fractions of A, T, G and C bases at each base location of reads or in the entire dataset.
# agtc_cont <- atgcContent(files_list)

# return the proportion of mapped reads included in a SAM/BAM file. 
# For paired end reads, it can return the proportion of mapped fragments (ie. read pairs).
# prop_mapped <- propmapped(files_list)

# Extract quality strings and convert them to Phred scores for generating boxplots.
# Quality scores give the probabilities of read bases being incorrectly called
# x <- qualityScores(filename, input_format = "gzFASTQ", offset = 33,nreads = 10000)

# Arguments
# filename = a character string giving the name of an input file containing sequence reads.
# input_format = a character string specifying format of the input file. gzFASTQ (gzipped FASTQ) by default. Acceptable formats include gzFASTQ, FASTQ, SAM and BAM. Character string is case insensitive.
# offset = a numeric value giving the offset added to the base-calling Phred scores. Possible values include 33 and 64. By default, 33 is used.
# nreads = a numeric value giving the number of reads from which quality scores are extracted. 10000 by default. A value of -1 indicates that quality scores will be extracted from all the reads.
