#!/usr/bin/env bash
#Cutadapt
################
cutadapt ./1_ORIGIN/SRR3182445_1.fastq -u -1 -o ./3_CUTADAPT/  >./log/SRR3182445_1_Cutadapt_standard_output 2>./log/_Cutadapt_log
cutadapt ./1_ORIGIN/SRR3182445_2.fastq -u -1 -o ./3_CUTADAPT/  >./log/SRR3182445_2_Cutadapt_standard_output 2>./log/_Cutadapt_log

#BWA&GATK
################Read Group ID
    normal_sample_ID=
    reference_genome=
    gatk_reference=
    known_indel=
    mills=
    dbsnp=

bwa mem -M -R '@RG\tID:${read_group_ID}\tSM:NORMAL\tLB:WXS\tPL:Illumina' ${reference_genome} ./3_CUTADAPT/${normal_sample_ID}_1_trimmed.fastq ./3_CUTADAPT/${normal_sample_ID}_2_trimmed.fastq  | samtools sort -o ./4_BWA_GATK/${normal_sample_ID}.bam -
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

###################################Read Group ID
    tumor_sample_ID=

bwa mem -M -R '@RG\tID:${read_group_ID}\tSM:TUMOR\tLB:WXS\tPL:Illumina' ${reference_genome} ./3_CUTADAPT/${tumor_sample_ID}_1_trimmed.fastq ./3_CUTADAPT/${tumor_sample_ID}_2_trimmed.fastq  | samtools sort -o ./4_BWA_GATK/${tumor_sample_ID}.bam -
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

#Mutation Calling
###################
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
--normalBam=./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--tumorBam=./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--referenceFasta=${gatk_reference} \
--runDir=${tumor_sample_ID}_manta

${tumor_sample_ID}_manta/runWorkflow.py -m local

configureStrelkaSomaticWorkflow.py --exome \
--normalBam=./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--tumorBam=./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--referenceFasta=${gatk_reference} \
--runDir=${tumor_sample_ID}_strelka \
--indelCandidates=${tumor_sample_ID}_manta/results/variants/candidateSmallIndels.vcf.gz

${tumor_sample_ID}_strelka/runWorkflow.py -m local

cat ${tumor_sample_ID}_strelka/results/variants/somatic.indels.vcf.gz|gzip -dc > ./5_VARIANT_CALLING/${tumor_sample_ID}_ms.indels.vcf
cat ${tumor_sample_ID}_strelka/results/variants/somatic.snvs.vcf.gz|gzip -dc > ./5_VARIANT_CALLING/${tumor_sample_ID}_ms.snvs.vcf
tar -f - -cv ${tumor_sample_ID}_manta | xz -T0 -9 >${tumor_sample_ID}_manta.txz
tar -f - -cv ${tumor_sample_ID}_strelka | xz -T0 -9 >${tumor_sample_ID}_strelka.txz
rm -rf ${tumor_sample_ID}_manta/ ${tumor_sample_ID}_strelka/

#Mutect
mutect --analysis_type MuTect \
--reference_sequence ${gatk_reference} \
--input_file:normal ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam \
--input_file:tumor ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam \
--vcf ${tumor_sample_ID}_mutect1.vcf \
--out ${tumor_sample_ID}_mutect1.out

xz -9 -T0 ${tumor_sample_ID}_mutect1.out
mv ${tumor_sample_ID}_mutect1.vcf ./5_VARIANT_CALLING/${tumor_sample_ID}_mutect1.vcf

#Mutect2
gatk Mutect2 -R ${gatk_reference} -I ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam -normal ${mutect2_normal_sample_ID_name} -I ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam -tumor ${mutect2_tumor_sample_ID_name} -O ${tumor_sample_ID}_mutect2_raw.vcf

gatk FilterMutectCalls -R ${gatk_reference} -V ${tumor_sample_ID}_mutect2_raw.vcf -O ${tumor_sample_ID}_mutect2.vcf

mv ${tumor_sample_ID}_mutect2.vcf ./5_VARIANT_CALLING/${tumor_sample_ID}_mutect2.vcf

#Shimmer
shimmer ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam --ref ${gatk_reference} --outdir ${tumor_sample_ID}_shimmer

perl add_readct_shimmer.pl ${tumor_sample_ID}_shimmer/som_counts.bh.txt ${tumor_sample_ID}_shimmer/somatic_diffs.vcf ${tumor_sample_ID}_shimmer/somatic_diffs_readct.vcf

cat ${tumor_sample_ID}_shimmer/somatic_diffs_readct.vcf|gawk '{OFS="\t";FS="\t";if ($0 ~ /^#/){print $0}else{if ($6>=25){$7="PASS"} print $0}}'>${tumor_sample_ID}_shimmer.vcf

mv ${tumor_sample_ID}_shimmer.vcf ./5_VARIANT_CALLING/${tumor_sample_ID}_shimmer.vcf

tar -f - -cv ${tumor_sample_ID}_shimmer |xz -9 -T0>${tumor_sample_ID}_shimmer.txz
rm -rf ${tumor_sample_ID}_shimmer/

#varscan
samtools mpileup -f ${gatk_reference} -q 1 -B ./4_BWA_GATK/${normal_sample_ID}_FINAL.bam ./4_BWA_GATK/${tumor_sample_ID}_FINAL.bam | \
varscan somatic Opt: ${tumor_sample_ID}_varscan-qbrc --output-vcf 1 --mpileup 1

varscan processSomatic ${tumor_sample_ID}_varscan-qbrc.indel.vcf --min-tumor_sample_ID-freq 0.01
varscan processSomatic ${tumor_sample_ID}_varscan-qbrc.snp.vcf --min-tumor_sample_ID-freq 0.01

mv ${tumor_sample_ID}_varscan-qbrc.snp.Somatic.hc.vcf ./5_VARIANT_CALLING/${tumor_sample_ID}_varscan-qbrc.snvs.vcf
mv ${tumor_sample_ID}_varscan-qbrc.indel.Somatic.hc.vcf ./5_VARIANT_CALLING/${tumor_sample_ID}_varscan-qbrc.indels.vcf

#Filtering
####################
    normal_sample_ID=
    tumor_sample_ID=
    Annovar=

Rscript BEFORE_ANNOVAR_FILTERING.R \
/gpfsdata/home/kangjiangn/Proj/Test6/Shi_et_al_2018_Retest \
./5_VARIANT_CALLING/${normal_sample_ID}_lofreq-call.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_lofreq-call.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_mutect1.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_mutect2.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_shimmer.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_ms.indels.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_ms.snvs.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_varscan-qbrc.indels.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_varscan-qbrc.snvs.vcf \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations.txt

table_annovar.pl ./6_FILTERING/${tumor_sample_ID}_somatic_mutations.txt  ${Annovar} \
-buildver hg38 \
-out ./6_FILTERING/${tumor_sample_ID}_somatic_mutations \
-remove -nastring . \
-protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all \
-operation g,f,f,f,f,f

Rscript AFTER_ANNOVAR_FILTERING.R \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations.txt \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations.hg38_multianno.txt \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations_after_R.txt

annotate_variation.pl -geneanno \
-dbtype refGene -buildver hg38 \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations_after_R.txt ${Annovar}

coding_change.pl --includesnp --alltranscript \
--newevf ./6_FILTERING/${tumor_sample_ID}_somatic_FINAL.txt \
./6_FILTERING/${tumor_sample_ID}_somatic_mutations_after_R.txt.exonic_variant_function \
${Annovar}/hg38_refGene.txt \
${Annovar}/hg38_refGeneMrna.fa

cat ./6_FILTERING/${tumor_sample_ID}_somatic_FINAL.txt|cut -f 4,5,7,8|sed 's/chr//'>"${tumor_sample_ID}_somatic_FINAL.loc"

mv ${tumor_sample_ID}_somatic_FINAL.loc ./6_FILTERING/${tumor_sample_ID}_somatic_FINAL.loc
