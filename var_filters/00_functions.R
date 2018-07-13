# -----------------------------------
# Function script for my Sepidermidis Project
# Author: Nick Brazeau
# -----------------------------------


# -----------------------------------
# Functions for 02_vcf_filter
# -----------------------------------

# vcf to tidy dataframe
vcf_info_field_to_tidy <-  function(vcf){
  infolist <- c("AF", "DP", "QD", "MQ", "SOR")
  infolist <- lapply(infolist, 
                     function(x){vcfR::extract_info_tidy(vcf, info_fields = x)})
  infodf <- plyr::join_all(infolist, by = "Key", type = "left")
  
  if(typeof(infodf$AF) == "character"){
    infodf$AF <- as.numeric(infodf$AF) # odd default in extract_info -- but character to numeric is fine in R
  }
  return(infodf)
}



#-----------------------------------------------------
# Read VCF and Go to VCF of Filtered Sites
#------------------------------------------------------
vcffilter <- function(vcffile = NULL, 
                      formatGQ=30, 
                      infoMQ=50,
                      infoQD=25,
                      infoSOR=2,
                      infoAF = 0.05,
                      infoDP = NULL, 
                      prop.loci.missing = 0.05,
                      biallelic = TRUE,
                      SNPs = FALSE){
  
  require(vcfR)
  require(tidyverse)
  
  # -----------------------------------------------------
  # Read and check input
  #------------------------------------------------------
  vcf <- vcfR::read.vcfR(file=vcffile)
  if(biallelic==TRUE){
    vcf <- vcf[vcfR::is.biallelic(vcf)] # subset to biallelic
  }
  if(SNPs==TRUE){
    vcf <-vcfR::extract.indels(vcf, return.indels = F) # subset to SNPs
  }
  # store loci objects on info fields 
  infolist <- c("AF", "DP", "QD", "MQ", "SOR")
  infolist <- lapply(infolist, 
                     function(x){vcfR::extract_info_tidy(vcf, info_fields = x)})
  infodf <- plyr::join_all(infolist, by = "Key", type = "left")
  
  if(typeof(infodf$AF) == "character"){
    infodf$AF_t <- as.numeric(infodf$AF) # odd default in extract_info -- but character to numeric is fine in R
  }
  
  #--------------------------------------------------------
  # filter loci
  #--------------------------------------------------------
  if(!is.null(infoAF)){
    infodf <- infodf %>% 
      dplyr::filter(AF >= infoAF & AF <= 1-infoAF)
  } 
  if(!is.null(infoDP)){
    DPpercentile <- quantile(infodf$DP, c(infoDP, 1-infoDP))
    infodf <- infodf %>% 
      dplyr::filter(DP >= DPpercentile[1] & DP <= DPpercentile[2])
  }  
  if(!is.null(infoMQ)){
    infodf <- infodf %>% 
      dplyr::filter(MQ >= infoMQ)
  }  
  if(!is.null(infoQD)){
    infodf <- infodf %>% 
      dplyr::filter(QD >= infoQD)
  }  
  if(!is.null(infoSOR)){
    infodf <- infodf %>% 
      dplyr::filter(SOR <= infoSOR)
  }  
  passedloci <- infodf$Key
  
  
  #--------------------------------------------------------
  # filter sample-level GQ
  #--------------------------------------------------------
  # store format and filter fields
  if(!is.null(formatGQ)){
    vcfsample <- vcfR::extract.gt(vcf, element =  "GQ", as.numeric = T)
    vcfsample[vcfsample < formatGQ] <- NA
    vcfsample <- cbind.data.frame(vcf@gt[,1], vcfsample) # need to add format column so matrix vcf@gt will match

  #--------------------------------------------------------
  # Subset by loci, GQ
  #--------------------------------------------------------
  vcf@gt[is.na(vcfsample)] <- NA
  
  locimissingness <- apply(vcf@gt, 1, function(x){sum(is.na(x))})
  locimissingness <- which(locimissingness == (ncol(vcf@gt)-1)) # in case a whole loci excluded by GQ
  passedloci <- !passedloci %in% locimissingness
  
  # filter bad loci
  vcf@gt <- vcf@gt[passedloci,]
  
  #--------------------------------------------------------
  # Drop samples with prop of loci missing
  #--------------------------------------------------------
  sample.prop.loci.missing <- colSums(is.na(vcf@gt))/nrow(vcf@gt)
  vcf@gt <- vcf@gt[, c(prop.loci.missing > sample.prop.loci.missing)]
  } else{
    
    # filter bad loci (if GQ not included)
    vcf@gt <- vcf@gt[passedloci,]
  }
  
  fix <- as.matrix(vcfR::getFIX(vcf, getINFO = T)[passedloci,])
  gt <- as.matrix(vcf@gt)
  meta <- append(vcf@meta, "##Additional Filters provided by polyIBD filter tools")
  meta <- append(vcf@meta, paste("Some samples may have been filtered by polyIBD filter tools. The new sample count is:", ncol(gt)-1))
  
  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  
  newvcfR
  
}


#-----------------------------------------------------
# Read vcfR object to segregated sites
#------------------------------------------------------

vcfR2segsites <- function(vcfRobj = NULL, err = 0.025){
  f_vcfAD <- vcfR::AD_frequency(vcfR::extract.gt(vcfRobj, element="AD"))
  if(any(is.na(f_vcfAD))){
    stop("There are NA values in your vcf matrix. Do you have sites that are not biallelic? If so, parsing will not work.")
  }
  segsites <- apply(f_vcfAD, 1, function(x){
  !all(x <= (err + median(x)) & (x >= (err - median(x))))
      }
  )
  
  vcfRobj@gt <- vcfRobj@gt[segsites,]

  fix <- as.matrix(vcfR::getFIX(vcfRobj, getINFO = T)[segsites,])
  gt <- as.matrix(vcfRobj@gt)
  meta <- append(vcfRobj@meta, paste("##Additional Filters for segregating sites, such that within-sample AF must vary more than:", err))

  # Setting class based off of vcfR documentation https://github.com/knausb/vcfR/blob/master/R/AllClass.R
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  
  newvcfR
  
}
