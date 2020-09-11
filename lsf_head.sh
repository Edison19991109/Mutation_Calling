#!/usr/bin/env bash
# BSUB -R "span[ptile=1]"
# BSUB -o
# BSUB -n 1
set -ue
WDIR=
cd "${WDIR}"
line=
thisdate="$(date +%Y-%m-%d_%H-%M)"
. "${HOME}"/lib/libdo # You should Change this line to path of `libdo` of your YuZJLab LinuxMiniPrograms installation.
LIBDO_LOG="${WDIR}/${line}_${thisdate}.log"
LIBDO_LOG_MODE=4
# You should modify following lines:
HG="${HOME}/References/GATK/hg38"
ref="${HG}/Homo_sapiens_assembly38.fasta"
onekg_snp="${HG}/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
mills="${HG}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
dbsnp="${HG}/dbsnp_146.hg38.vcf.gz"
omni="${HG}/1000G_omni2.5.hg38.vcf.gz"
hapmap="${HG}/hapmap_3.3.hg38.vcf.gz"
k_indel=" -known ${onekg_snp} -known ${mills} -L ${HG}/wgs_calling_regions.hg38.interval_list "
ANNOVAR_DB="${HOME}/References/Annovar/humandb"
bwaref="${HOME}/References/GATK_hg38/hg38_bwa_ref/hg38.fa"
