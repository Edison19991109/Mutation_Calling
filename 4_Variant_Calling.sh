#!/usr/bin/env bash
# BSUB -R "span[ptile=1]"
# BSUB -o ./log/Cutadapt_log
# BSUB -e ./log/Cutadapt_error
# BSUB -n 1
    gatk_reference=
    tumor_sample_ID=
    normal_sample_ID=
    mutect2_normal_sample_ID_name=
    mutect2_tumor_sample_ID_name=

#lofreq-call
lofreq call -C 7 -s --sig 0.1 --bonf 1 -f ${gatk_reference} --call-indels -o ./5_VARIANT_CALLING/${tumor_sample_ID}_lofreq-call.vcf ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam
lofreq call -C 7 -s --sig 0.99 --bonf 1 -f ${gatk_reference} --call-indels -o ./5_VARIANT_CALLING/${normal_sample_ID}_lofreq-call.vcf ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam

#MS
configManta.py --exome \
--normal_sample_IDBam=./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--tumor_sample_IDBam=./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--referenceFasta=${gatk_reference} \
--runDir=${normal_sample_ID}_manta

${normal_sample_ID}_manta/runWorkflow.py -m local

configureStrelkaSomaticWorkflow.py --exome \
--normal_sample_IDBam=./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--tumor_sample_IDBam=./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--referenceFasta=${gatk_reference} \
--runDir=${normal_sample_ID}_strelka \
--indelCandidates=${normal_sample_ID}_manta/results/variants/candidateSmallIndels.vcf.gz

${normal_sample_ID}_strelka/runWorkflow.py -m local

cat ${normal_sample_ID}_strelka/results/variants/somatic.indels.vcf.gz|gzip -dc > ./5_VARIANT_CALLING/${normal_sample_ID}_ms.indels.vcf
cat ${normal_sample_ID}_strelka/results/variants/somatic.snvs.vcf.gz|gzip -dc > ./5_VARIANT_CALLING/${normal_sample_ID}_ms.snvs.vcf
tar -f - -cv ${normal_sample_ID}_manta | xz -T0 -9 >${normal_sample_ID}_manta.txz
tar -f - -cv ${normal_sample_ID}_strelka | xz -T0 -9 >${normal_sample_ID}_strelka.txz
rm -rf ${normal_sample_ID}_manta/ ${normal_sample_ID}_strelka/

#Mutect
mutect --analysis_type MuTect \
--reference_sequence ${gatk_reference} \
--input_file:normal_sample_ID ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--input_file:tumor_sample_ID ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--vcf ${normal_sample_ID}_mutect1.vcf \
--out ${normal_sample_ID}_mutect1.out

xz -9 -T0 ${normal_sample_ID}_mutect1.out
mv ${normal_sample_ID}_mutect1.vcf ./5_VARIANT_CALLING/${normal_sample_ID}_mutect1.vcf

#Mutect2
gatk Mutect2 -R ${gatk_reference} -I ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam -normal_sample_ID ${mutect2_normal_sample_ID_name} -I ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam -tumor_sample_ID ${mutect2_tumor_sample_ID_name} -O ${normal_sample_ID}_mutect2_raw.vcf

gatk FilterMutectCalls -R ${gatk_reference} -V ${normal_sample_ID}_mutect2_raw.vcf -O ${normal_sample_ID}_mutect2.vcf

mv ${normal_sample_ID}_mutect2.vcf ./5_VARIANT_CALLING/${normal_sample_ID}_mutect2.vcf

#Shimmer
shimmer ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam --ref ${gatk_reference} --outdir ${normal_sample_ID}_shimmer

perl add_readct_shimmer.pl ${normal_sample_ID}_shimmer/som_counts.bh.txt ${normal_sample_ID}_shimmer/somatic_diffs.vcf ${normal_sample_ID}_shimmer/somatic_diffs_readct.vcf

cat ${normal_sample_ID}_shimmer/somatic_diffs_readct.vcf|gawk '{OFS="\t";FS="\t";if ($0 ~ /^#/){print $0}else{if ($6>=25){$7="PASS"} print $0}}'>VCF/${normal_sample_ID}_shimmer.vcf

mv ${normal_sample_ID}_shimmer.vcf ./5_VARIANT_CALLING/${normal_sample_ID}_shimmer.vcf

tar -f - -cv ${normal_sample_ID}_shimmer |xz -9 -T0>${normal_sample_ID}_shimmer.txz
rm -rf ${normal_sample_ID}_shimmer/

#varscan
samtools mpileup -f ${gatk_reference} -q 1 -B ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \|\
varscan somatic Opt: ${normal_sample_ID}_varscan-qbrc --output-vcf 1 --mpileup 1

varscan processSomatic ${normal_sample_ID}_varscan-qbrc.indel.vcf --min-tumor_sample_ID-freq 0.01
varscan processSomatic ${normal_sample_ID}_varscan-qbrc.snp.vcf --min-tumor_sample_ID-freq 0.01

mv ${normal_sample_ID}_varscan-qbrc.snp.Somatic.hc.vcf ./5_VARIANT_CALLING/${normal_sample_ID}_varscan-qbrc.snvs.vcf
mv ${normal_sample_ID}_varscan-qbrc.indel.Somatic.hc.vcf ./5_VARIANT_CALLING/${normal_sample_ID}_varscan-qbrc.indels.vcf
