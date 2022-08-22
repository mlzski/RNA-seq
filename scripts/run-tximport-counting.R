# alternative version of 'run_DEG_swish.R' to be used only for counts generation
# it extracts counts directly from quants.sf files, either on transcript- and gene-level (using tximport package)

# alternatively, it could be also done using 'tximeta' package to get output as SummarizedExperiment object

if (!requireNamespace("tximport", quietly = TRUE)) BiocManager::install("tximport")
if (!requireNamespace("tximeta", quietly = TRUE)) BiocManager::install("tximeta")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
suppressMessages(library(tximport))
suppressMessages(library(tximeta))
suppressMessages(library(SummarizedExperiment))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

args <- commandArgs(trailingOnly = TRUE)

# for testing ONLY
# args <- c("/nobackup/ummz/analysis-May-22/read-quantification/output-quantif/coh-1",
#         "/nobackup/ummz/analysis-May-22/counts/coh-1",
#         "/nobackup/ummz/analysis-May-22/counts/samples-info-coh-1.txt") 

# args <- c("/nobackup/ummz/analysis-May-22/read-quantification/output-quantif/both", "/nobackup/ummz/analysis-May-22/counts-both-lvl", "/nobackup/ummz/analysis-May-22/read-quantification/output-quantif/samples-both.txt") 

# NOTE: both-levels means starting from transcript-level and then summarizing to gene-level

if (length(args)!=3) {
         stop("ERROR: 3 arguments must be supplied:
                 \n(1) INPUT directory,
                 \n(2) OUTPUT directory,
                 \n(3) sample information file", call.=FALSE)
}

# define whether output files should be saved or not [TRUE / FALSE]
output_save <- TRUE

# define directory with data (INPUT)
dir_in <- args[1]

# define directory for results (OUTPUT)
dir_out <- args[2]

# load table data
coldata <- read.table(args[3])

# get list of all files and add it to coldata table
coldata$files <- file.path(dir_in, coldata$V1, "quant.sf")
all(file.exists(coldata$files))

########################## to be removed ##########################

# -------------- tximeta --------------

# load quantification data 
se <- tximeta(coldata)

# get counts matrix
cts_tximeta <- assays(se)[["counts"]]

# save counts matrix as .csv
write.csv(cts_tximeta, file = file.path(args[2], "counts_transcript-level-tximeta.csv"))
cat("Created: ",  file.path(args[2], "counts_transcript-level-tximeta.csv"), "\n")
  
# save SummarizedExperiment object from tximeta
saveRDS(cts_tximeta, file = file.path(args[2], "obj_tximport_transcript-level-tximeta.rds"))
cat("Created: ", file.path(args[2], "obj_tximport_transcript-level-tximeta.rds"), "\n")

####################################################################

### transcript-level ###
  
# import quantification data
obj_tximport <- tximport(coldata$files, type = "salmon", txOut = TRUE)

# get counts matrix
cts_tximport <- obj_tximport$counts

# add column names
colnames(cts_tximport) <- coldata$IDs

# save counts matrix as .csv
write.csv(cts_tximport, file = file.path(args[2], "counts_transcript-level.csv"))
cat("Created: ", file.path(args[2], "counts_transcript-level.csv"), "\n")

# save list object from tximport
saveRDS(obj_tximport, file = file.path(args[2], "obj_tximport_transcript-level.rds"))
cat("Created: ", file.path(args[2], "obj_tximport_transcript-level.rds"), "\n")

### gene-level summarization ###

# make a data.frame called tx2gene with two columns: transcriptID and geneID (column order matters)
# NOTE: transcriptID must be the same one used in the abundance files

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

# -------------------- GET SUMMARY TABLE -------------------- #
# use the addIds function from tximeta to add gene symbols, By specifying gene=TRUE, 
# this will use the gene ID to match to gene symbols for all of the transcripts.
se <- addIds(se, "SYMBOL", gene=TRUE)

# save the summary table
cts_transcript_summary <- rowRanges(se)
write.table(as.data.frame(cts_transcript_summary), file = file.path(args[2], "transcript_summary.csv"), sep=";")

cat("Created: ", file.path(args[2], "transcript_summary.csv"), "\n")

# gene-level

# summarize all of the data to the gene level
gse <- summarizeToGene(se)

# get counts matrix
cts_gene <- assays(gse)[["counts"]]

# save table with results 
write.csv(cts_gene, file = file.path(args[2], "genes_counts.csv"))

gse <- addIds(gse, "SYMBOL", gene=TRUE)

# save the summary table
cts_gene_summary <- rowRanges(gse)
write.table(as.data.frame(cts_gene_summary), file = file.path(args[2], "gene_summary.csv"), sep=";")

cat("Finished ! ! !")