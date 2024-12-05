#!/bin/bash

# File path and parameters
MIN_READ_LENGTH=20
N_THREAD= # add number of threads
FASTQ_PATH= # add path to folder of FASTQ files
TRIMMED_FASTQ_PATH= # add path to folder where to generate trimmed FASTQ files
LOG_FILES_PATH= # add path to folder of where to generate log files

# Trimming and adapter removing
for FASTQ in $(ls $FASTQ_PATH | grep '_R1.fastq.gz$'); do
    SAMPLE_NAME=$(basename $FASTQ '_R1.fastq.gz')
    READ1="$FASTQ_PATH/${SAMPLE_NAME}_R1.fastq.gz"
    READ2="$FASTQ_PATH/${SAMPLE_NAME}_R2.fastq.gz"
    echo "Trimming $SAMPLE_NAME" &> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    ~/Desktop/TrimGalore-0.6.6/trim_galore --paired --dont_gzip --no_report_file --length $MIN_READ_LENGTH -j $N_THREAD -o $TRIMMED_FASTQ_PATH/ $READ1 $READ2 &>> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
done

