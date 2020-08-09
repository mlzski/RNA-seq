### this script runs Picard Tools to locate and tag duplicate reads in a BAM or SAM file

if [ $# != 3 ] ; then
    echo -e "ERROR: 3 arguments are required: \
	(1) path to data folder [BAM files], \
	(2) path to output directory, \
	(3) SGE_TASK_ID argument for array jobs ... Exiting"
    exit 1
fi    

# assign variables
data_dir=$1   #                                                   | 11026_S12_Aligned.sortedByCoord.out.bam
out_dir=$2    # /nobackup/ummz/analyses/rerun_Ian/rerun_1/pic     | 11026_S12_Aligned.sortedByCoord.out.bam
#METRICS_FILE # /nobackup/ummz/analyses/rerun_Ian/rerun_1/pic     | 11026_S12_Aligned.sortedByCoord.out.bam.metrics.txt

#core_name=`echo $data_dir | rev | cut -d'/' -f 1 | rev`

bamFile=$(ls $data_dir/*.bam* | sed -n -e "$SGE_TASK_ID p")

core_name=$(ls $data_dir/*.bam* | rev | cut -d'/' -f 1 | rev | sed -n -e "$SGE_TASK_ID p")

java -jar /home/home02/ummz/tools/picard/build/libs/picard.jar MarkDuplicates \
	INPUT=$bamFile \					# => (String) 			One or more input SAM or BAM files to analyze. Must be coordinate sorted.
	OUTPUT=$out_dir/$core_name \				# => (File)			The output file to write marked records to Required.
	METRICS_FILE=$out_dir/${core_name}.metrics.txt \	# => (File) 			File to write duplication metrics to Required.
	TAGGING_POLICY=All \					# => (DuplicateTaggingPolicy)	Determines how duplicate types are recorded in the DT optional attribute. instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). 
	TMP_DIR=$out_dir/temp \					# => (File)			Default value: null. This option may be specified 0 or more times.
	VALIDATION_STRINGENCY=LENIENT 				# => (ValidationStringency)	Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded. Default value: STRICT. This option can be set to 'null' to clear the default value. Possible values: {STRICT, LENIENT, SILENT}

###############################################################################################################################################################
# parameters guide: (copied from Picard documentation)

# default options that were set in Ian's file:
#	MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=50000 \		# => (Integer)		This option is obsolete. ReadEnds will always be spilled to disk. Default value: 50000. This option can be set to 'null' to clear the default value.
#	MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \		# => (Integer)		Maximum number of file handles to keep open when spilling read ends to disk. Set this number a little lower than the per-process maximum number of file that may be open. This number can be found by executing the 'ulimit -n' command on a Unix system. Default value: 8000. This option can be set to 'null' to clear the default value.
#	SORTING_COLLECTION_SIZE_RATIO=0.25 \			# => (Double)		This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some of the sorting collections. If you are running out of memory, try reducing this number. Default value: 0.25. This option can be set to 'null' to clear the default value.
#	ASSUME_SORTED=false \					# => (Boolean)		If true, assume that the input file is coordinate sorted even if the header says otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false} Cannot be used in conjuction with option(s) ASSUME_SORT_ORDER (ASO)
#	DUPLICATE_SCORING_STRATEGY=SUM_OF_BASE_QUALITIES \	# => (ScoringStrategy)	The scoring strategy for choosing the non-duplicate among candidates. Default value: SUM_OF_BASE_QUALITIES. This option can be set to 'null' to clear the default value. Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}
# 	PROGRAM_RECORD_ID=MarkDuplicates \			# => (String)		The program record ID for the @PG record(s) created by this program. Set to null to disable PG record creation. This string may have a suffix appended to avoid collision with other program record IDs. Default value: MarkDuplicates. This option can be set to 'null' to clear the default value.
#	PROGRAM_GROUP_NAME=MarkDuplicates \			# => (String)		Value of PN tag of PG record to be created. Default value: MarkDuplicates. This option can be set to 'null' to clear the default value.
#	OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \			# => (Integer)		The maximum offset between two duplicate clusters in order to consider them optical duplicates. The default is appropriate for unpatterned versions of the Illumina platform. For the patterned flowcell models, 2500 is moreappropriate. For other platforms and models, users should experiment to find what works best. Default value: 100. This option can be set to 'null' to clear the default value.
#	REMOVE_SEQUENCING_DUPLICATES=false \			# => (Boolean)		If true remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process instead of the library preparation process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and this option is ignored. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
#	REMOVE_DUPLICATES=false \				# => (Boolean)		If true do not write duplicates to the output file instead of writing them with appropriate flags set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
#	VERBOSITY=INFO \					# => (LogLevel)		Control verbosity of logging. Default value: INFO. This option can be set to 'null' to clear the default value. Possible values: {ERROR, WARNING, INFO, DEBUG}
#	QUIET=false \						# => (Boolean)		Whether to suppress job-summary info on System.err. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# 	COMPRESSION_LEVEL=5 \					# => (Integer) 		Compression level for all compressed files created (e.g. BAM and GELI). Default value: 5. This option can be set to 'null' to clear the default value.
#	MAX_RECORDS_IN_RAM=500000 \				# => (Integer)		When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed. Default value: 500000. This option can be set to 'null' to clear the default value.
#	CREATE_INDEX=false \					# => (Boolean)		Whether to create a BAM index when writing a coordinate-sorted BAM file. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
#	CREATE_MD5_FILE=false \					# => (Boolean)		Whether to create an MD5 digest for any BAM or FASTQ files created. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# 	GA4GH_CLIENT_SECRETS=client_secrets.json		# => (String)		Google Genomics API client_secrets.json file path. Default value: client_secrets.json. This option can be set to 'null' to clear the default value.

# other options (not used above)
# BARCODE_TAG 			(String)	Barcode SAM tag (ex. BC for 10X Genomics) Default value: null.
# READ_ONE_BARCODE_TAG 		(String)	Read one barcode SAM tag (ex. BX for 10X Genomics) Default value: null.
# READ_TWO_BARCODE_TAG 		(String)	Read two barcode SAM tag (ex. BX for 10X Genomics) Default value: null.
# TAG_DUPLICATE_SET_MEMBERS 	(Boolean)	If a read appears in a duplicate set, add two tags. The first tag, DUPLICATE_SET_SIZE_TAG (DS), indicates the size of the duplicate set. The smallest possible DS value is 2 which occurs when two reads map to the same portion of the reference only one of which is marked as duplicate. The second tag, DUPLICATE_SET_INDEX_TAG (DI), represents a unique identifier for the duplicate set to which the record belongs. This identifier is the index-in-file of the representative read that was selected out of the duplicate set. Default value: false. This option can be set to 'null' to clear the default value. Possible values: {true, false}
# ASSUME_SORT_ORDER 		(SortOrder)	If not null, assume that the input file has this order even if the header says otherwise. Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate, unknown} Cannot be used in conjuction with option(s) ASSUME_SORTED (AS)
# PROGRAM_GROUP_VERSION 	(String)	Value of VN tag of PG record to be created. If not specified, the version will be detected automatically. Default value: null.
# PROGRAM_GROUP_COMMAND_LINE 	(String)	Value of CL tag of PG record to be created. If not supplied the command line will be detected automatically. Default value: null.
# COMMENT 			(String)	Comment(s) to include in the output file's header. Default value: null. This option may be specified 0 or more times.
# READ_NAME_REGEX 		(String)	Regular expression that can be used to parse read names in the incoming SAM file. Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate. These values are used to estimate the rate of optical duplication in order to give a more accurate estimated library size. Set this option to null to disable optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate. The regular expression should contain three capture groups for the three variables, in order. It must match the entire read name. Note that if the default regex is specified, a regex match is not actually done, but instead the read name is split on colon character. For 5 element names, the 3rd, 4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8), the 5th, 6th, and 7th elements are assumed to be tile, x and y values. Default value: . This option can be set to 'null' to clear the default value.

# NOTICE: the field 'ESTIMATED_LIBRARY_SIZE' will only be filled for PE data (empty for SE)
# "The estimated number of unique molecules in the library based on PE duplication."

# tool used: 'MarkDuplicates' => identifies duplicate reads. This tool locates and tags duplicate reads in a BAM or SAM file.
# where duplicate reads are defined as originating from a single fragment of DNA. 
# Duplicates can arise during sample preparation e.g. library construction using PCR. 
# Duplicates can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.

# The 'MarkDuplicates' tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. 
# An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. 
# After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).

# INPUT: (either coordinate-sorted or query-sorted inputs)
# => when the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. 
# => when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.

# OUTPUT:
# => a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. (Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.)
# => a metrics file indicating the numbers of duplicates for both single- and paired-end reads.

# usage example:
# java -jar picard.jar MarkDuplicates \
#      	I=input.bam \			# => bam file obtained from STAR
#      	O=marked_duplicates.bam \	# => new bam file generated by Picard
#      	M=marked_dup_metrics.txt	# => file with metrics
