rm(list = ls())
require(maftools)
options(stringsAsFactors = F)
## annovar
annovar.laml <- annovarToMaf(annovar = "test.hg38_multianno.txt", 
                             refBuild = 'hg38',
                             tsbCol = 'Tumor_Sample_Barcode', 
                             table = 'refGene',
                             MAFobj = T)
write.table(annovar.laml,file='all_mutations',quote=F,sep="\t",row.names = F)
## gatk
# gatk.laml = read.maf(maf = 'gatk/gatk4.1.4.0_merge.maf')

# library(data.table)
# tmp=fread('./7.annotation/funcatator/funcatator_merge.maf')
# gatk.laml = read.maf(maf = tmp)
# 
# ## vep
# vep.laml = read.maf(maf = './7.annotation/vep/vep_merge.maf')
# ## for vep.laml
# library(stringr)
# vep.laml@data$Protein_Change = paste0("p.",
#                                       str_sub(vep.laml@data$Amino_acids,1,1),
#                                       vep.laml@data$Protein_position,
#                                       str_sub(vep.laml@data$Amino_acids,3,3))

## save Rdata
save(annovar.laml, file = 'laml.Rdata')

## Summary
laml=annovar.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)

laml=gatk.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)

laml=vep.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)

rm(list = ls())
require(maftools)
options(stringsAsFactors = F)
load(file = 'laml.Rdata')
#laml=c(annovar.laml,gatk.laml,vep.laml)

## mafsummary
#anno=c('annovar','gatk','vep')
  #i = 1
  pdf('plotmafSummary.pdf',width = 12,height = 8)
  plotmafSummary(maf = annovar.laml,
                 rmOutlier = TRUE,
                 showBarcodes = T,
                 textSize = 0.4,
                 addStat = 'median',
                 dashboard = TRUE,
                 titvRaw = FALSE)
  dev.off()
#Oncoplot
  pdf('Oncoplottop30.pdf',width=10,height=10)
  oncoplot(maf = annovar.laml,
           top = 30,
           fontSize = 0.5,
           sampleOrder = annovar.laml@clinical.data$Tumor_Sample_Barcode,
           showTumorSampleBarcodes = T)
  dev.off()
  