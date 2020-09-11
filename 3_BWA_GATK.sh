#!/usr/bin/env bash
# BSUB -R "span[ptile=1]"
# BSUB -o ./log/Cutadapt_log
# BSUB -e ./log/Cutadapt_error
# BSUB -n 1
    normal_sample_ID=
    reference_genome=
    gatk_reference=
    known_indel=
    mills=
    dbsnp=
    read_group_ID=

bwa mem -M -R '@RG'\""\tID:${read_group_ID}\tSM:NORMAL\tLB:WXS\tPL:Illumina"\" ${reference_genome} ./3_CUTADAPT/${normal_sample_ID}_1.fastq ./3_CUTADAPT/${normal_sample_ID}_2.fastq  \| samtools sort -o ./4_BWA_GATK/${normal_sample_ID}.bam -
samtools index ./4_BWA_GATK/${normal_sample_ID}.bam

picard MarkDuplicates \
    INPUT=./4_BWA_GATK/${normal_sample_ID}.bam \
    OUTPUT=./4_BWA_GATK/${normal_sample_ID}_dupmark.bam \
    REMOVE_SEQUENCING_DUPLICATES=true \
    METRICS_FILE=./4_BWA_GATK/${normal_sample_ID}_dupmark_metrics.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=STRICT

gatk3 -T RealignerTargetCreator \
-R ${gatk_reference} ${known_indel} \
-o ./4_BWA_GATK/${normal_sample_ID}_intervals.list \
-I ./4_BWA_GATK/${normal_sample_ID}_dupmark.bam

gatk3 -T IndelRealigner \
--filter_bases_not_stored \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R ${gatk_reference} ${known_indel} \
-targetIntervals ./4_BWA_GATK/${normal_sample_ID}_intervals.list \
-I ./4_BWA_GATK/${normal_sample_ID}_dupmark.bam \
-o ./4_BWA_GATK/${normal_sample_ID}_realigned.bam

gatk BaseRecalibrator \
-R ${gatk_reference} \
-known-sites ${mills} \
-known-sites ${dbsnp} \
-I ./4_BWA_GATK/${normal_sample_ID}_realigned.bam -O ./4_BWA_GATK/${normal_sample_ID}_bqsr

gatk ApplyBQSR \
-R ${gatk_reference} \
-I ./4_BWA_GATK/${normal_sample_ID}_realigned.bam \
-O ./4_BWA_GATK/${normal_sample_ID}_bqsr.bam \
-bqsr ./4_BWA_GATK/${normal_sample_ID}_bqsr

samtools view ./4_BWA_GATK/${normal_sample_ID}_bqsr.bam -q 30 -b -o ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam

samtools index ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam

###################################
    tumor_sample_ID=
    read_group_ID=

bwa mem -M -R '@RG'\""\tID:${read_group_ID}\tSM:TUMOR\tLB:WXS\tPL:Illumina"\" ${reference_genome} ./3_CUTADAPT/${tumor_sample_ID}_1.fastq ./3_CUTADAPT/${tumor_sample_ID}_2.fastq  \| samtools sort -o ./4_BWA_GATK/${tumor_sample_ID}.bam -
samtools index ./4_BWA_GATK/${tumor_sample_ID}.bam

picard MarkDuplicates \
    INPUT=./4_BWA_GATK/${tumor_sample_ID}.bam \
    OUTPUT=./4_BWA_GATK/${tumor_sample_ID}_dupmark.bam \
    REMOVE_SEQUENCING_DUPLICATES=true \
    METRICS_FILE=./4_BWA_GATK/${tumor_sample_ID}_dupmark_metrics.txt \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=STRICT

gatk3 -T RealignerTargetCreator \
-R ${gatk_reference} ${known_indel} \
-o ./4_BWA_GATK/${tumor_sample_ID}_intervals.list \
-I ./4_BWA_GATK/${tumor_sample_ID}_dupmark.bam

gatk3 -T IndelRealigner \
--filter_bases_not_stored \
--disable_auto_index_creation_and_locking_when_reading_rods \
-R ${gatk_reference} ${known_indel} \
-targetIntervals ./4_BWA_GATK/${tumor_sample_ID}_intervals.list \
-I ./4_BWA_GATK/${tumor_sample_ID}_dupmark.bam \
-o ./4_BWA_GATK/${tumor_sample_ID}_realigned.bam

gatk BaseRecalibrator \
-R ${gatk_reference} \
-known-sites ${mills} \
-known-sites ${dbsnp} \
-I ./4_BWA_GATK/${tumor_sample_ID}_realigned.bam -O ./4_BWA_GATK/${tumor_sample_ID}_bqsr

gatk ApplyBQSR \
-R ${gatk_reference} \
-I ./4_BWA_GATK/${tumor_sample_ID}_realigned.bam \
-O ./4_BWA_GATK/${tumor_sample_ID}_bqsr.bam \
-bqsr ./4_BWA_GATK/${tumor_sample_ID}_bqsr

samtools view ./4_BWA_GATK/${tumor_sample_ID}_bqsr.bam -q 30 -b -o ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam

samtools index ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam