### this script runs featureCounts to count reads from STAR output or Picard output(name_Aligned.sortedByCoord.out.bam) 
# all files (samples) are launched in parallel and produces one output (.csv) file per each sample 

# install (if necessary) and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager",repos = "http://cran.us.r-project.org")
if (!requireNamespace("Rsubread", quietly = TRUE)) BiocManager::install("Rsubread")
library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)

# for running in R-studio
# args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam", 
#          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/5_/postprocessed")

if (length(args)!=5) {
  stop("5 arguments must be supplied: 
  	\n(1) running option [dups or nodups],
	\n(2) running mode parameter [SE or PE],
	\n(3 - input) path to directory with data, 
	\n(4 - annotation) path to GTF annotation files and
	\n(5 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[3], sep="\n")

cat("Directory for results (OUT): ")
cat(args[5], sep="\n")

setwd(args[5])

# retrive built-in annotation file (other genomes available: mm9, mm10, hg19 and hg38 (NCBI RefSeq))
# ann <- getInBuiltAnnotation(annotation = "hg38")
# These In-built annotations were downloaded from NCBI RefSeq database and then 
# => adapted by merging overlapping exons from the same gene to form a set of disjoint exons for each gene. 
# => Genes with the same Entrez gene identifiers were also merged into one gene.
# => Each row in the annotation represents an exon of a gene. 
# => There are five columns in the annotation data including Entrez gene identifier (GeneID), chromosomal name (Chr), chromosomal start position(Start), chromosomal end position (End) and strand (Strand).

# define path to GTF annotation files
ann <- args[4]

# create a list with file names to be processed:
files_list <- list.files(args[3], pattern = ".bam$", full.names = TRUE)

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

# run featureCounts() command [SE - single end | PE - paired end] AND
# write output tables with counts and statistics
# URL to parameters: https://www.rdocumentation.org/packages/Rsubread/versions/1.22.2/topics/featureCounts 

# INPUT:
# (1) input a set of files that contain read mapping results
# (2) an annotation file that includes genomic features. (either in SAM and BAM format, which is detected automatically) 
# NOTICE: Input reads can be name-sorted or location-sorted. 

fc_res <- featureCounts(files_list, 
                        annot.ext=ann,
			isGTFAnnotationFile = TRUE, 	
                        allowMultiOverlap=TRUE, 	
                        fraction=TRUE, 	
                        nthreads=8, 
                        ignoreDup=parDups, 		
                        isPairedEnd=parMode, 		
                        useMetaFeatures=TRUE)		

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

