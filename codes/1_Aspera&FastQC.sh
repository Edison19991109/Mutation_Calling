#!/usr/bin/env bash
    sample_header=
    sample_id=

    for i in ${sample_id};
    do
    var=`expr ${i} % 10`
    ascp -P 33001 -QT -l 300m  -i asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR181/00${var}/${sample_header}${i}/${sample_header}${i}.fastq.gz ./1_ORIGIN >./log/${sample_header}${i}_ASCP_standard_output 2>./log/${sample_header}${i}_RAW_log
    gunzip ./1_ORIGIN/${sample_header}${i}.fastq.gz 2>./log/${sample_header}${i}_gunzip_log
    done;

    for i in ${sample_id};
    do
    fastqc --threads 3 -o ./2_FASTQC ./1_ORIGIN/${sample_header}${i}.fastq >./log/${sample_header}${i}_FASTQC_standard_output 2>./log/${sample_header}${i}_FASTQC_log;
    done;