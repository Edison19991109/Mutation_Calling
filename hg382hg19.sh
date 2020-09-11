#!/usr/bin/env bash
cat ../all.conf | grep -v germline | cut -d " " -f 1 | while read id; do
    CrossMap.py vcf "${HOME}"/References/CHAIN/hg38ToHg19.over.chain.gz ${id}_varscan.indel.vcf "${HOME}"/References/GATK/hg19/ucsc.hg19.fasta ${id}_varscan_FINAL_converted.vcf
done
