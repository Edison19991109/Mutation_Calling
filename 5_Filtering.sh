#!/usr/bin/env bash
    normal_sample_ID=
    tumor_sample_ID=
    Annovar=

DO Rscript exec/BEFORE_ANNOVAR_FILTERING.R \
./5_VARIANT_CALLING/${normal_sample_ID}_lofreq-call.vcf \
./5_VARIANT_CALLING/${tumor_sample_ID}_lofreq-call.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_mutect1.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_mutect2.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_shimmer.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_ms.indels.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_ms.snvs.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_varscan-qbrc.indels.vcf \
./5_VARIANT_CALLING/${normal_sample_ID}_varscan-qbrc.snvs.vcf \
./6_FILTERING/${normal_sample_ID}_somatic_mutations.txt

DO table_annovar.pl ./6_FILTERING/${normal_sample_ID}_somatic_mutations.txt  ${Annovar} \
-buildver hg38 \
-out ./6_FILTERING/${normal_sample_ID}_somatic_mutations \
-remove -nastring . \
-protocol refGene,ljb26_all,cosmic70,esp6500siv2_all,exac03,1000g2015aug_all \
-operation g,f,f,f,f,f

DO Rscript exec/AFTER_ANNOVAR_FILTERING.R \
./6_FILTERING/${normal_sample_ID}_somatic_mutations.txt \
./6_FILTERING/${normal_sample_ID}_somatic_mutations.hg38_multianno.txt \
./6_FILTERING/${normal_sample_ID}_somatic_mutations_after_R.txt

DO annotate_variation.pl -geneanno \
-dbtype refGene -buildver hg38 \
./6_FILTERING/${normal_sample_ID}_somatic_mutations_after_R.txt ${Annovar}

DO coding_change.pl --includesnp --alltranscript \
--newevf ./6_FILTERING/${normal_sample_ID}_somatic_FINAL.txt \
./6_FILTERING/${normal_sample_ID}_somatic_mutations_after_R.txt.exonic_variant_function \
./6_FILTERING/${Annovar}/hg38_refGene.txt \
./6_FILTERING/${Annovar}/hg38_refGeneMrna.fa

cat ./6_FILTERING/${normal_sample_ID}_somatic_FINAL.txt|cut -f 4,5,7,8|sed 's/chr//'>loc/"${normal_sample_ID}_somatic_FINAL.loc"

mv ${normal_sample_ID}_somatic_FINAL.loc ./6_FILTERING/${normal_sample_ID}_somatic_FINAL.loc
