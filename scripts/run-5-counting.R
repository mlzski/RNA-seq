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

fc_res <- featureCounts(files_list, 			# 
                        annot.ext=ann, 			# a character string giving name of a user-provided annotation file, one of these files: hg38.ensGene.gtf, hg38.knownGene.gtf, hg38.ncbiRefSeq.gtf (recommended), hg38.refGene.gtf
			isGTFAnnotationFile = TRUE, 	# logical indicating whether the annotation provided via the annot.ext argument is in GTF format or not. FALSE by default.
                        allowMultiOverlap=TRUE, 	# logical indicating if a read is allowed to be assigned to more than one feature (or meta-feature) if it is found to overlap with more than one feature (or meta-feature). FALSE by default.
                        fraction=TRUE, 			# logical indicating if fractional counts are produced for multi-mapping reads and/or multi-overlapping reads. FALSE by default.
                        nthreads=8, 			#
                        ignoreDup=parDups, 		# logical indicating whether reads marked as duplicates should be ignored. FALSE by default. Read duplicates are identified using bit Ox400 in the FLAG field in SAM/BAM files. The whole fragment (read pair) will be ignored if paired end.
                        isPairedEnd=parMode, 		# logical indicating if counting should be performed on read pairs or reads. FALSE by default. If TRUE, read pairs will be counted instead of individual reads.
                        useMetaFeatures=TRUE)		# logical indicating whether the read summarization should be performed at the feature level (eg. exons) or meta-feature level (eg. genes). If TRUE, features in the annotation (each row is a feature) will be grouped into meta-features, using the GTF.attrType attribute in the GTF-format annotation file, and then reads will be assiged to the meta-features instead of the features.

# parameters set to default:
# 			countMultiMappingReads=TRUE, 	# logical indicating if multi-mapping reads/fragments should be counted, TRUE by default. ‘NH’ tag is used to located multi-mapping reads in the input BAM/SAM files.                        

# parameters that I'm not sure about
#			strandSpecific=2		# an integer vector indicating if strand-specific read counting should be performed. Length of the vector should be either 1 (meaning that the value is applied to all input files), or equal to the total number of input files provided. Each vector element should have one of the following three values: 0 (unstranded), 1 (stranded) and 2 (reversely stranded). Default value of this parameter is 0 (ie. unstranded read counting is performed for all input files).

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
