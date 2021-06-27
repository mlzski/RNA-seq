
# alternative version of 'run_DEG_swish.R' to be used only for counts generation, on transcript- and gene-level
# in includes both methods: tximport() and txmeta()

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")
if (!requireNamespace("fishpond", quietly = TRUE)) BiocManager::install("fishpond")
if (!requireNamespace("tximeta", quietly = TRUE)) BiocManager::install("tximeta")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("tximport", quietly = TRUE)) BiocManager::install("tximport")

suppressMessages(library(SummarizedExperiment))
suppressMessages(library(fishpond))
suppressMessages(library(tximeta))
suppressMessages(library(tximport))
suppressMessages(library(org.Hs.eg.db))

# get working directory to recognise the machine
w_dir <- getwd()

# create a shortcut for the OneDrive directory where all files are stored
if(startsWith(w_dir, "/Users/michal")){           
  main_dir <- "/Users/michal/Documents/OneDrive - University of Leeds"    # on my mac
} else if (startsWith(w_dir, "/Users/ummz")) {    
  main_dir <- "/Users/ummz/Documents/OneDrive - University of Leeds"      # on uni mac    
} else if (startsWith(w_dir, "/nobackup/ummz")) {
  main_dir <- "/nobackup/ummz"						  # on ARC4
} else {
  print("Unrecognised machine.")
}

args <- commandArgs(trailingOnly = TRUE)

# for testing ONLY
#args <- c("/nobackup/ummz/analyses/run_16_Apr21/quants_1-22",
#         "/nobackup/ummz/swish_v2/output_1-22",
#         "/nobackup/ummz/swish_v2/samples.txt") 

#args <- c("/nobackup/ummz/analyses/run_15_Mar21/quants_all", 
#	  "/nobackup/ummz/swish_v3/test",
# 	  "/nobackup/ummz/swish_v3/samples_all.txt")
 
if (length(args)!=3) {
         stop("ERROR: 3 arguments must be supplied: 
                 \n(1) INPUT directory,
                 \n(2) OUTPUT directory,
                 \n(3) annotation file", call.=FALSE)
}

# define wheather output files should be saved or not [TRUE / FALSE]
output_save <- TRUE

# define directory with data (INPUT)
dir_in <- args[1]

# define directory for results (OUTPUT)
dir_out <- args[2]

# load table data
coldata <- read.table(args[3], header = TRUE)

# define the nem to be used and create a new folder for results of the current run
# colnames(coldata)
run_feat <- "counts"
  
# check if directory for this feature already exists
if(dir.exists(file.path(dir_out, run_feat)) == FALSE){
  # if not, create it
  dir.create(file.path(dir_out, run_feat))
  paste0("Created new directory: ", file.path(dir_out, run_feat))
} else {
  cat("Output directory already exists")
}

# get list of all files and add it to coldata table
coldata$files <- file.path(dir_in, coldata$names, "quant.sf")
all(file.exists(coldata$files))

# load in the quantification data with tximeta
se <- tximeta(coldata)

# transcript-level
y <- se

# get counts matrix
cts_transcript_tximeta <- assays(y)[["counts"]]

# save transcript-level matrix as .csv
write.csv(cts_transcript_tximeta, file = file.path(dir_out, run_feat, paste0("transcript_level_tximeta_", run_feat, ".csv")))

# import quantification data with tximport
txi.tx <- tximport(coldata$files, type = "salmon", txOut = TRUE)

# get counts matrix
cts_transcript_tximport <- txi.tx$counts

colnames(cts_transcript_tximport) <- coldata$names

# save transcript-level matrix as .csv
write.csv(cts_transcript_tximport, file = file.path(dir_out, run_feat, paste0("transcript_level_tximport_", run_feat, ".csv")))
 
# use the addIds function from tximeta to add gene symbols. 
# By specifying gene=TRUE, this will use the gene ID to match to gene symbols for all of the transcripts.
y <- addIds(y, "SYMBOL", gene=TRUE)

# save the summary table
cts_transcript_summary <- rowRanges(y)
write.table(as.data.frame(cts_transcript_summary), file = file.path(dir_out, run_feat, "transcript_level_summary.csv"), sep=";")

# gene-level

# summarize all of the data to the gene level
gse <- summarizeToGene(se)

gy <- gse

# get counts matrix
cts_gene <- assays(gy)[["counts"]]

# save table with results 
write.csv(cts_gene, file = file.path(dir_out, run_feat, paste0("gene_level_", run_feat, ".csv")))

gy <- addIds(gy, "SYMBOL", gene=TRUE)

# save the summary table
cts_gene_summary <- rowRanges(gy)
write.table(as.data.frame(cts_gene_summary), file = file.path(dir_out, run_feat, "gene_level_summary.csv"), sep=";")

cat("Finished ! ! !")

