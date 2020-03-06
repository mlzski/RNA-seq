###################################################################################################### 
# this script runs featureCounts to count reads from STAR output (name_Aligned.sortedByCoord.out.bam)
###################################################################################################### 

# INPUT:
# (1) input a set of files that contain read mapping results
# (2) an annotation file that includes genomic features. (either in SAM and BAM format, which is detected automatically) 

# Input reads can be name-sorted or location-sorted. 
# Users do not need to resort the reads before feeding them to featureCounts.

# use getInBuiltAnnotation() to retrived built-in annotation files: for genomes mm9, mm10, hg19 and hg38 (NCBI RefSeq)

# install (if necessary) and load packages


library(Rsubread)

args <- commandArgs(trailingOnly = TRUE)

# for running in R-studio
# args <- c("/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/bam", 
#          "/Users/ummz/OneDrive - University of Leeds/ANALYSES/results_run_I_Nov19/4_alignement/postprocessed")

if (length(args)!=2) {
  stop("2 arguments must be supplied: \n(1 - input) path to directory with data and \n(2 - output) path where output files should be stored", call.=FALSE)
}

cat("Directories with data (IN): ")
cat(args[1], sep="\n")

cat("Directory for results (OUT): ")
cat(args[2], sep="\n")

setwd(args[1])

# retrive built-in annotation file
ann <- getInBuiltAnnotation(annotation = "hg38")

# create a list with file names to be processed:
files_list <-
  
  
fc_SE <- featureCounts(files_list, annot.ext=ann, nthreads = 1)

# for paired-end (the bam file needs to ogrinate from paired-end mode of alignment)
# fc_PE <- featureCounts("alignResultsPE.BAM", annot.ext=ann, isPairedEnd=TRUE)

# OUTPUTS:
# fc_SE$counts      => a very long vector with counts for all genes
# fc_SE$annotation  => a table with all genes (rows) 6 columns for: "GeneID", "Chr", "Start", "End", "Strand", "Length"
# fc_SE$targets     => file name(s) used
# fc_SE$stat        => a table with statistics (14 metrics)

# concatenate all fc_SE$counts


# Extract quality strings and convert them to Phred scores for generating boxplots.
# Quality scores give the probabilities of read bases being incorrectly called
x <- qualityScores(filename, input_format = "gzFASTQ", offset = 33,nreads = 10000)

# Arguments
# filename = a character string giving the name of an input file containing sequence reads.
# input_format = a character string specifying format of the input file. gzFASTQ (gzipped FASTQ) by default. Acceptable formats include gzFASTQ, FASTQ, SAM and BAM. Character string is case insensitive.
# offset = a numeric value giving the offset added to the base-calling Phred scores. Possible values include 33 and 64. By default, 33 is used.
# nreads = a numeric value giving the number of reads from which quality scores are extracted. 10000 by default. A value of -1 indicates that quality scores will be extracted from all the reads.


# return fractions of A, T, G and C bases at each base location of reads or in the entire dataset.
# atgcContent("8546_S3_L005_Aligned.sortedByCoord.out.bam")

# return the proportion of mapped reads included in a SAM/BAM file. For paired end reads, it can return the proportion of mapped fragments (ie. read pairs).
# propmapped("8546_S3_L005_Aligned.sortedByCoord.out.bam")






