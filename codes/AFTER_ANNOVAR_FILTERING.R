#Read in the original file before ANNOVAR(the file name should look like mutation.txt).
options(scipen=999)
args <- commandArgs(trailingOnly = TRUE)
mutations <- tryCatch(read.table(file = args[1],stringsAsFactors = F,sep="\t",colClasses=c(rep(NA,3),rep("character",2),rep(NA,5))),
  error=function(e) {
    as.data.frame(matrix(0,ncol=10,nrow=0))
  }
)

#Assign the column name of mutation data matrix.
colnames(mutations) <- c("Chr","Start","End","Ref","Alt","Caller","Normal_ref","Normal_alt",
                      "Tumor_ref","Tumor_alt")

#Read in the file name that have got the ANNOVAR annotation (the file name should look like multianno.txt).
##Fill in the file name.
annotations <- read.table(file = args[2],stringsAsFactors = F,sep="\t",header=T)
annotations <- cbind(mutations,annotations[,-c(1:5)])

as_numeric<-function(x)
{
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
}

# human/PDX sample
if ("cosmic70" %in% colnames(annotations))
{
  suppressWarnings({annotations$cosmic70 <- as_numeric(annotations$cosmic70)})
  suppressWarnings({annotations$esp6500siv2_all <- as_numeric(annotations$esp6500siv2_all)})
  suppressWarnings({annotations$ExAC_ALL <- as_numeric(annotations$ExAC_ALL)})
  suppressWarnings({annotations$X1000g2015aug_all <- as_numeric(annotations$X1000g2015aug_all)})
  ##suppressWarnings({annotations$EAS.sites.2015_08 <- as_numeric(annotations$EAS.sites.2015_08)})
  suppressWarnings({annotations$ExAC_EAS <- as_numeric(annotations$ExAC_EAS)})
  
    annotations <- annotations[annotations$esp6500siv2_all<0.01,]
    annotations <- annotations[annotations$ExAC_ALL<0.01,]
    annotations <- annotations[annotations$X1000g2015aug_all<0.01,]
    ##annotations <- annotations[annotations$ExAC_EAS<0.01,]
    ##annotations <- annotations[annotations$EAS.sites.2015_08<0.01,]
  
  effect_field <- c("SIFT_pred","Polyphen2_HVAR_pred","cosmic70","esp6500siv2_all",
                  "ExAC_ALL","X1000g2015aug_all")
}

annotations <- annotations[,c("Chr","Start","End","Ref","Alt","Caller","Normal_ref","Normal_alt",
 "Tumor_ref","Tumor_alt","Func.refGene","Gene.refGene","ExonicFunc.refGene",
 "AAChange.refGene",effect_field)]

#Input the output file name.
write.table(annotations, file= args[3], sep= "\t", row.names=F, quote=F)
