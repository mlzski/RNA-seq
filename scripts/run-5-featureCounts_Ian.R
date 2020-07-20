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

if (length(args)!=5) {
  stop("5 arguments must be supplied: 
	\n(1) running option [dups or nodups],
	\n(2) running mode [SE or PE],
	\n(3 - input) path to directory with data, 
	\n(4 - annotation) path to GTF annotation files
	\n(5 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[3], sep="\n")

cat("Directory for results (OUT): ")
cat(args[5], sep="\n")

setwd(args[5])

# define path to GTF annotation files
ann <- args[4]

# create a list with file names to be processed:
files_list <- list.files(args[3], pattern = ".bam$", full.names = TRUE)

# run featureCounts() command [dups - with duplicates | nodups - without duplicates] AND
# write output tables with counts and statistics

# URL to parameters: https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts 

# loops for running modes
if(args[1] == "dups"){
	parDups <- as.logical(FALSE)
} else if (args[1] == "nodups"){
	parDups <- as.logical(TRUE)
} else {
	stop("ERROR: running option must be defined as either dups or nodups", call.=FALSE)	
}

if(args[2] == "SE"){
	parMode <- as.logical(FALSE)
} else if (args[2] == "PE"){
	parMode <- as.logical(TRUE)
} else {
	stop("ERROR: running mode must be defined as either SE or PE", call.=FALSE)	
}
	
fc_res <- featureCounts(files_list, 
			annot.ext=ann, 
			isGTFAnnotationFile = 
			TRUE, useMetaFeatures=TRUE, 
			allowMultiOverlap=TRUE, 
			countMultiMappingReads=TRUE, 
			fraction=TRUE, nthreads=8, 
			ignoreDup=parDups, 
			isPairedEnd=parMode, 
			strandSpecific=2)

if(args[1] == "dups"){
	write.csv(fc_res$counts, file="all_counts_dups.csv")
	write.csv(fc_res$stat, file="all_stats_dups.csv")
	write.csv(fc_res$annotation, file="all_annotations_dups.csv")

} else if (args[1] == "nodups"){
	write.csv(fc_res$counts, file="all_counts_nodups.csv")
	write.csv(fc_res$stat, file="all_stats_nodups.csv")
	write.csv(fc_res$annotation, file="all_annotations_nodups.csv")
	
} else {
	stop("ERROR: running option must be defined as either dups or nodups", call.=FALSE)	
}

proc.time()
