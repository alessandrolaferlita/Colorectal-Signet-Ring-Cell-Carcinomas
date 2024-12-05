#!/bin/bash

# File path and parameters
N_THREAD= # add number of thread
INDEXED_GENOME= # add path to indexed genome generate by HISAT2
BASENAME_INDEXED_GENOME= # add basename of the files of the indexed genome
TRIMMED_FASTQ_PATH= # add path to folder of trimmed fastq
SAM_FILES_PATH= # add path to folder of where to generate SAM files
BAM_FILES_PATH= # add path to folder of where to generate BAM files
SORTED_BAM_FILES_PATH= # add path to folder of where to generate sorted BAM files
LOG_FILES_PATH= # add path to folder of where to generate log files

# Alignment
for TRIMMED_FASTQ in $(ls $TRIMMED_FASTQ_PATH | grep '_R1_val_1.fq$'); do
    SAMPLE_NAME=$(basename $TRIMMED_FASTQ '_R1_val_1.fq')
    SAM_NAME=$SAMPLE_NAME".sam"
    READ1="$TRIMMED_FASTQ_PATH/${SAMPLE_NAME}_R1_val_1.fq"
    READ2="$TRIMMED_FASTQ_PATH/${SAMPLE_NAME}_R2_val_2.fq"
    echo "Aligning sample $SAMPLE_NAME" &> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    hisat2 -p $N_THREAD -x $INDEXED_GENOME/$BASENAME_INDEXED_GENOME -1 $READ1 -2 $READ2 -S $SAM_FILES_PATH/$SAM_NAME &>> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    # Conversion SAM to BAM
    BAM_NAME=$SAMPLE_NAME".bam"
    echo "Converting $SAM_NAME" &>> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    samtools view -bS $SAM_FILES_PATH/$SAM_NAME > $BAM_FILES_PATH/$BAM_NAME
    rm -R $SAM_FILES_PATH/*
    # Sorting BAM
    SORTED_BAM_NAME=$SAMPLE_NAME'_sorted.bam'
    echo "Sorting $SAMPLE_NAME" &>> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    samtools sort --threads $N_THREAD $BAM_FILES_PATH/$BAM_NAME -o $SORTED_BAM_FILES_PATH/$SORTED_BAM_NAME
    rm -R $BAM_FILES_PATH/*
    # Index BAM
    BAI_NAME=$SAMPLE_NAME'_sorted.bam.bai'
    echo "Indexing $SAMPLE_NAME" &>> $LOG_FILES_PATH/$SAMPLE_NAME'_log.txt'
    samtools index -@ $N_THREAD $SORTED_BAM_FILES_PATH/$SORTED_BAM_NAME $SORTED_BAM_FILES_PATH/$BAI_NAME
done

echo "Finish!"

