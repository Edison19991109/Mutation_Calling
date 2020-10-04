#!/usr/bin/env bash
    ${files_that_needed_to_be_converted}=
    ${files_that_are_converted}=
    CrossMap.py vcf /References/CHAIN/hg38ToHg19.over.chain.gz ${files_that_needed_to_be_converted} /References/GATK/hg19/ucsc.hg19.fasta ${files_that_are_converted}
