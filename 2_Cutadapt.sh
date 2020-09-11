#!/usr/bin/env bash
cutadapt ./2_FASTQC/ -u -u -a -o ./3_CUTADAPT/ _trimmed.fastq  >./log/ _Cutadapt_standard_output 2>./log/ _Cutadapt_log