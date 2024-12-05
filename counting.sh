#!/bin/bash

# File path and parameters
N_THREAD= # add number of threads
BAM_FILES_PATH= # add path to folder containing the sorted BAM files
RAW_COUNTS_PATH= # add path to folder where to generate the output of featureCounts
GTF_FILE= # add path to gencode.v43.annotation.gtf (this file can be download from GenCode website)

# Counting
featureCounts $BAM_FILES_PATH/*bam -a $GTF_FILE -p -g gene_id -T $N_THREAD -o $RAW_COUNTS_PATH/raw_counts_table.txt

