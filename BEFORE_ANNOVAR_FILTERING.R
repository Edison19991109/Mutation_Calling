#############################  read vcfs  #########################
##Fancy methods checking the installed libraries and load them.
rm(list = ls())
name_pkg <- c("statmod")
bool_nopkg <- !name_pkg %in% rownames(installed.packages())
if(sum(bool_nopkg) > 0){
  install.packages(name_pkg[bool_nopkg])
}
invisible(lapply(name_pkg, library, character.only = T)) # load multiple packages


##This works under the Linux bash background.
##If you run it under the Rstudio, then you can just skip it.
options(scipen=999)
args <- commandArgs(trailingOnly = TRUE)

#Input the working direcotry.
wd <- args[1]

#Input the name of the each mutation calling program file.
#Example input: lofreq_n.vcf
lofreq_n_ReadIn <- args[2]
#Example input: lofreq_t.vcf
lofreq_t_ReadIn <- args[3]

##Read in the vcf file
mutect1_ReadIn <- args[4]

mutect2_ReadIn <- args[5]

#Example input: somatic_diffs.readct.vcf
shimmer_ReadIn <- args[6]

#Example input: passed.somatic.indels.vcf
strelka_indel_ReadIn <- args[7]
#Example input: passed.somatic.snvs.vcf
strelka_snp_ReadIn <- args[8]

#Example input: varscan.indel.Somatic.hc.vcf
varscan_indel_ReadIn <- args[9]
#Example input: varscan.snp.Somatic.hc.vcf
varscan_snp_ReadIn <- args[10]

out_file <- args[11]

#Set working directory.
setwd(wd)


##This works under the Rstuido background
##If you run it under the Linux bash background, then you can just skip it.

##Set up the working directory.
# working_directory <- c()
# setwd(working_directory)

#########################  Functions Set Up  ########################
##This part should be run first.

##Read in the vcf.
read_vcf <- function(file)
{
  cat(paste("Reading",file,"\n"))
  ##If error happens during read in file, return an empty vcf file.
  x <- tryCatch({
        read.table(file,stringsAsFactors = F)
      }, error = function(e) {
         
        ## Give an empty vcf file
        cat(paste("Warning: failed to read",file,
                  ". Maybe this caller didn't find any variants\n"))
        data.frame(V1=character(0),V2=numeric(0),V4=character(0),V5=character(0),V7=character(0),
                  V8=character(0),normal_ref=numeric(0),normal_alt=numeric(0),
                  tumor_ref=numeric(0),tumor_alt=numeric(0),stringsAsFactors = F)
      
      })
  x
}

##The extract_count function is used in the filter_vcf function.
extract_count <- function(fields,info,field)
{
  if (length(info)==0) {return(numeric(0))}
  split_info <- strsplit(info,":")
  as.numeric(sapply(1:length(info),
                    function(i) sub(",.*","",split_info[[i]][fields[[i]]==field])))
}

##Filtering the vcf as criteria defined.
##Only use the somatic filtering option in this function.
filter_vcf <- function(vcf,caller,type="somatic")
{
  ##Used for debugging.
  #browser()

  cat(paste("Filtering",caller,type,"\n"))
  if (dim(vcf)[1]==0) {return(vcf)}
  
  ##Doubt remaining for change the content of V8
  ##May need to rewrite.
  vcf$V8 <- caller
  
  ##Doubt remain for the useage at least in Mutect.
  vcf <- vcf[!grepl(",",vcf$V5),]
  
  ##Filter and leave those passing the mutation calling program filtering criteria.
  vcf <- vcf[vcf$V7=="PASS",]
  
  ##Calculate the stopping position.
  vcf$V3 <- vcf$V2+nchar(vcf$V4)-1
  
  ##Normal and mutation hit is not null.
  vcf <- vcf[vcf$V10!="." & vcf$V11!=".",]
  
  ##Read in and record the relative content of V10 as well as V11.
  fields <- strsplit(vcf$V9,":")
  
  # extract read count
  ##Extract depends on caller.
  if (caller=="mutect")
  {
    ##Return a null data frame if there is no mutation left.
    if (dim(vcf)[1]==0)
    {
      vcf$normal_ref=vcf$normal_alt=vcf$tumor_ref=vcf$tumor_alt=numeric(0)
      return(vcf)
    }
    
    normal_ct <- strsplit(sapply(strsplit(vcf$V11,":"),function(x) x[2]),",")
    vcf$normal_ref <- as.numeric(sapply(normal_ct,function(x) x[1]))
    vcf$normal_alt <- as.numeric(sapply(normal_ct,function(x) x[2]))
    
    tumor_ct <- strsplit(sapply(strsplit(vcf$V10,":"),function(x) x[2]),",")
    vcf$tumor_ref <- as.numeric(sapply(tumor_ct,function(x) x[1]))
    vcf$tumor_alt <- as.numeric(sapply(tumor_ct,function(x) x[2]))
  }else if (caller %in% c("speedseq","shimmer","varscan","strelka"))
  {
    if (caller %in% c("speedseq","shimmer"))
    {
      RO="RO";AO="AO"
    }else if (caller %in% c("varscan"))
    {
      RO="RD";AO="AD"
    }else if (caller %in% c("strelka"))
    {
      RO="TAR";AO="TIR"
    }
    
    vcf$normal_ref <- extract_count(fields,vcf$V10,RO)
    vcf$normal_alt <- extract_count(fields,vcf$V10,AO)
    vcf$tumor_ref <- extract_count(fields,vcf$V11,RO)
    vcf$tumor_alt <- extract_count(fields,vcf$V11,AO)
  }else if (caller=="lofreq") # lofreq, only for unmatched samples
  {
    tmp <- sapply(fields,function(x) strsplit(x,";")[[1]],simplify=F)
    total_count <- sapply(tmp,function(x) as.numeric(sub("DP=","",grep("DP=",x,value=T))))
    var_count <- round(total_count*
                         sapply(tmp,function(x) as.numeric(sub("AF=","",grep("AF=",x,value=T)))))
    vcf$normal_alt=vcf$tumor_alt=var_count
    vcf$normal_ref=vcf$tumor_ref=total_count-var_count
  }else if (caller=="strelka_snp") 
  {
    if (dim(vcf)[1]>0)
    {
      vcf$normal_alt=vcf$tumor_alt=vcf$normal_ref=vcf$tumor_ref=NA
      for (i in 1:dim(vcf)[1])
      {
        tmp <- strsplit(vcf$V10[i],":")[[1]]
        vcf$normal_ref[i] <- as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V4[i],"U",sep="")]))
        vcf$normal_alt[i] <- as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V5[i],"U",sep="")]))
        tmp <- strsplit(vcf$V11[i],":")[[1]]
        vcf$tumor_ref[i] <- as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V4[i],"U",sep="")]))
        vcf$tumor_alt[i] <- as.numeric(sub(".*,","",tmp[fields[[i]]==paste(vcf$V5[i],"U",sep="")]))
      }
    }else
    {
      vcf$normal_ref=vcf$normal_alt=vcf$tumor_ref=vcf$tumor_alt=numeric(0)
    }
  }else if (caller=="strelka_germline")
  {
    split_info <- strsplit(vcf$V10,":")
    tmp <- t(sapply(1:length(split_info),
                 function(i) as.numeric(strsplit(split_info[[i]][fields[[i]]=="AD"],",")[[1]])))
    vcf$normal_alt <- vcf$tumor_alt=tmp[,2]
    vcf$normal_ref <- vcf$tumor_ref=tmp[,1]
  }
  
  # Filter by read count and allele frequency
  ##The sum of normal read count and mutation read count greater than 5.
  ##Tumor mutation read count greater than 3.
  if (caller!="strelka_germline") 
  {
    vcf <- vcf[vcf$normal_ref+vcf$normal_alt>=5,]
    vcf <- vcf[vcf$tumor_ref+vcf$tumor_alt>=5,]
    vcf <- vcf[vcf$tumor_alt>=3,]
    if (type=="somatic")
    {
      vcf <- vcf[vcf$normal_alt/(vcf$normal_ref+vcf$normal_alt)<
                vcf$tumor_alt/(vcf$tumor_ref+vcf$tumor_alt)/2,]
      #vcf <- vcf[vcf$normal_alt/(vcf$normal_ref+vcf$normal_alt)<0.05,]
    }else
    {
      vcf <- vcf[vcf$normal_alt>=3,]
    } 
  }else # for tumor-only calling, make the calling super sensitive
  {
    vcf <- vcf[vcf$normal_ref+vcf$normal_alt>=3,]
    vcf <- vcf[vcf$normal_alt>=1,]
  }
  
  vcf
}




########################### Main Running Commands ##################
##Read in the mutect vcf file.
mutect1 <- read_vcf(mutect1_ReadIn)

##Filter the vcf file.
mutect1 <- filter_vcf(mutect1,"mutect")

##Read in the mutect2 vcf file.
mutect2 <- read_vcf(mutect2_ReadIn)
mutect2 <- filter_vcf(mutect2,"mutect")

#Read in the shimmer vcf file.
shimmer <- read_vcf(shimmer_ReadIn)
#It may need to change the V6 in the future.
shimmer <- shimmer[shimmer$V7=="PASS",]
if (dim(shimmer)[1]>0) {shimmer$V7="PASS"}
shimmer <- filter_vcf(shimmer,"shimmer")

#Read in the varscan file.
varscan_indel <- read_vcf(varscan_indel_ReadIn)
varscan_indel <- filter_vcf(varscan_indel,"varscan")
varscan_snp <- read_vcf(varscan_snp_ReadIn)
varscan_snp <- filter_vcf(varscan_snp,"varscan")
varscan <- rbind(varscan_indel,varscan_snp)

#Read in the strelka file
strelka_snp <- read_vcf(strelka_snp_ReadIn)
strelka_snp <- filter_vcf(strelka_snp,"strelka_snp")
if (dim(strelka_snp)[1]>0) {strelka_snp$V8="strelka"}
strelka_indel <- read_vcf(strelka_indel_ReadIn)
if (dim(strelka_indel)[1]>0) {strelka_indel$V7="PASS"}
strelka_indel <- filter_vcf(strelka_indel,"strelka")
strelka <- rbind(strelka_snp,strelka_indel)

#Read in the lofreq file.
lofreq_t <- read_vcf(lofreq_t_ReadIn)
lofreq_t$V9=lofreq_t$V10=lofreq_t$V11=lofreq_t$V8
lofreq_t <- filter_vcf(lofreq_t,"lofreq","germline")
lofreq_t$mutation <- paste(lofreq_t$V1,lofreq_t$V2,lofreq_t$V4,lofreq_t$V5)

lofreq_n <- read_vcf(lofreq_n_ReadIn)
lofreq_n$V9=lofreq_n$V10=lofreq_n$V11=lofreq_n$V8
lofreq_n <- filter_vcf(lofreq_n,"lofreq","germline")
lofreq_n$mutation <- paste(lofreq_n$V1,lofreq_n$V2,lofreq_n$V4,lofreq_n$V5)

lofreq_t$normal_ref <- lofreq_n$normal_ref[match(lofreq_t$mutation,lofreq_n$mutation)]
lofreq_t$normal_alt <- lofreq_n$normal_alt[match(lofreq_t$mutation,lofreq_n$mutation)]
lofreq <- lofreq_t[,colnames(lofreq_t)!="mutation"]
lofreq_germline <- lofreq # for germline filtering
lofreq <- lofreq[is.na(lofreq$normal_alt) | (lofreq$normal_alt/(lofreq$normal_ref+lofreq$normal_alt)<
                lofreq$tumor_alt/(lofreq$tumor_ref+lofreq$tumor_alt)/2),]
# lofreq <- lofreq[is.na(lofreq$normal_alt) | 
#                 (lofreq$normal_alt/(lofreq$normal_ref+lofreq$normal_alt)<0.05),]
lofreq <- lofreq[lofreq$V6>=10,]


##Combine the two replicate result and left with those mutations that both samples have.
vcf <- rbind(mutect1, mutect2, shimmer, varscan, strelka, lofreq)[,c("V1","V2","V4","V5","V8","normal_ref",
                                 "normal_alt","tumor_ref","tumor_alt")]
colnames(vcf)[1:5] <- c("chr","pos","ref","alt","caller")
vcf$variant <- paste(vcf$chr,vcf$pos,vcf$ref,vcf$alt)
if (dim(vcf)[1]==0) {stop("No valid mutations left in VCF file!\n")}

tmp1 <- aggregate(vcf$caller,by=list(vcf$variant),function(x) paste(x,collapse=","))
tmp2 <- aggregate(vcf[,c("normal_ref","normal_alt","tumor_ref","tumor_alt")],
                  by=list(vcf$variant),function(x) mean(x,na.rm=T))
vcf <- cbind(tmp1,tmp2[,-1])
vcf$chr <- sapply(strsplit(vcf$Group.1," "),function(x) x[1])
vcf$pos <- as.numeric(sapply(strsplit(vcf$Group.1," "),function(x) x[2]))
vcf$ref <- sapply(strsplit(vcf$Group.1," "),function(x) x[3])
vcf$alt <- sapply(strsplit(vcf$Group.1," "),function(x) x[4])
vcf$pos2 <- vcf$pos+nchar(vcf$ref)-1

##Important! A variant must have been found by >=3 caller
##In this script, I adjust the command to choose those both appear in the original sample and replicate,
##hence it left those appear with a comma in the middle.
vcf_left <- vcf[grepl(",",vcf$x),] 

##Prepare to write the file.
vcf_left <- vcf_left[,c("chr","pos","pos2","ref","alt","x","normal_ref","normal_alt","tumor_ref","tumor_alt")]
write.table(vcf_left,file=out_file,col.names = F,row.names = F,sep="\t",quote=F)


# ##Just for Mutect1.
# mutect_vcf <- mutect1[,c("V1","V2","V4","V5","V8","normal_ref",
#                          "normal_alt","tumor_ref","tumor_alt")]
# colnames(mutect_vcf)[1:5] <- c("chr","pos","ref","alt","caller")
# mutect_vcf$variant <- paste(mutect_vcf$chr,mutect_vcf$pos,mutect_vcf$ref,mutect_vcf$alt)
# if (dim(mutect_vcf)[1]==0) {stop("No valid mutations left in VCF file!\n")}
# mutect_vcf$pos2 <- mutect_vcf$pos+nchar(mutect_vcf$ref)-1
# mutect_vcf <- mutect_vcf[,c("chr","pos","pos2","ref","alt","normal_ref","normal_alt","tumor_ref","tumor_alt")]
# write.table(mutect_vcf,file="mutect_somatic_mutations.txt",col.names = F,row.names = F,sep="\t",quote=F)

# #Extract the vcf file of each program.
# extract_each_program <- function(vcf, program_name){
#   vcf_vcf <- vcf[,c("V1","V2","V4","V5","V8","normal_ref",
#                            "normal_alt","tumor_ref","tumor_alt")]
#   colnames(vcf_vcf)[1:5] <- c("chr","pos","ref","alt","caller")
#   vcf_vcf$variant <- paste(vcf_vcf$chr,vcf_vcf$pos,vcf_vcf$ref,vcf_vcf$alt)
#   if (dim(vcf_vcf)[1]==0) {stop("No valid mutations left in VCF file!\n")}
#   vcf_vcf$pos2 <- vcf_vcf$pos+nchar(vcf_vcf$ref)-1
#   vcf_vcf <- vcf_vcf[,c("chr","pos","pos2","ref","alt","normal_ref","normal_alt","tumor_ref","tumor_alt")]
#   filename <- paste(program_name, "_somatic_mutations.txt")
#   write.table(vcf_vcf,file=filename,col.names = F,row.names = F,sep="\t",quote=F)
# }

# extract_each_program(mutect1, 'mutect1')
# extract_each_program(mutect2, 'mutect2')
# extract_each_program(shimmer, 'shimmer')
# extract_each_program(varscan, 'varscan')
# extract_each_program(strelka, 'MS')
# extract_each_program(lofreq, 'lofreq')