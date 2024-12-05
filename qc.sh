#!/bin/bash

# File path and parameters
N_THREAD= # add number of threads
FASTQ_PATH= # add path to folder of FASTQ files
OUTPUT_PATH= # add path to folder of where to generate FASTQC reports

# QC
for FASTQ in $(ls $FASTQ_PATH/*); do
    SAMPLE_NAME=$(basename $FASTQ '.fq')
    echo "Analyzing sample $SAMPLE_NAME"
    fastqc -t $N_THREAD $FASTQ --outdir=$OUTPUT_PATH
done

echo "Finish!"

