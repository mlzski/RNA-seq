###################################################################################################### 
# this script runs featureCounts to count reads from STAR output (name_Aligned.sortedByCoord.out.bam)
# adapted from Ian
###################################################################################################### 

# INPUT:
# (1) input a set of files that contain read mapping results (bam files)
# (2) an annotation file that includes genomic features. (either in SAM and BAM format, which is detected automatically) 
# NOTICE: Input reads can be name-sorted or location-sorted. 

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!requireNamespace("Rsubread", quietly = TRUE)) BiocManager::install("Rsubread")
library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)!=3) {
  stop("3 arguments must be supplied: 
	\n(1) running mode parameter [dups or nodups], 
	\n(2 - input) path to directory with data, 
	\n(3 - annotation) path to GTF annotation files
	\n(4 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[2], sep="\n")

cat("Directory for results (OUT): ")
cat(args[4], sep="\n")

setwd(args[4])

# define path to GTF annotation files
ann <- args[3]

# create a list with file names to be processed:
files_list <- list.files(args[2], pattern = ".bam$", full.names = TRUE)

# run featureCounts() command [dups - with duplicates | nodups - without duplicates] AND
# write output tables with counts and statistics

if(args[1] == "dups"){
	fc_dups <- featureCounts(files_list, 
							 annot.ext=ann, 
							 isGTFAnnotationFile = TRUE, 
						   	 useMetaFeatures=TRUE, 
						   	 allowMultiOverlap=TRUE, 
						   	 countMultiMappingReads=TRUE,
						   	 fraction=TRUE, 
						   	 nthreads=8, 
						   	 ignoreDup=FALSE, 
						   	 isPairedEnd=TRUE,
						   	 strandSpecific=2)

	write.csv(fc_dups$counts, file="all_counts_dups.csv")
	write.csv(fc_dups$stat, file="all_stats_dups.csv")
	write.csv(fc_dups$annotation, file="all_annotations_dups.csv")

} else if (args[1] == "nodups"){
	fc_nodups <- featureCounts(files_list, 
							   annot.ext=ann, 
							   isGTFAnnotationFile = TRUE, 
						   	   useMetaFeatures=TRUE, 
						   	   allowMultiOverlap=TRUE, 
						   	   countMultiMappingReads=TRUE,
						   	   fraction=TRUE, 
						   	   nthreads=8, 
						   	   ignoreDup=TRUE, 
						   	   isPairedEnd=TRUE,
						   	   strandSpecific=2)

	write.csv(fc_nodups$counts, file="all_counts_nodups.csv")
	write.csv(fc_nodups$stat, file="all_stats_nodups.csv")
	write.csv(fc_nodups$annotation, file="all_annotations_nodups.csv")

} else {
	stop("ERROR: running mode parameter must be defined as either dups or nodups", call.=FALSE)	
}

proc.time()
