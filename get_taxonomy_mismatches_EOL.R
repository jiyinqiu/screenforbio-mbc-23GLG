#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

message("Welcome to get_taxonomy_mismatches.R")
message("")
message("Step1: Check for and load library(taxize)")
pkgLoad <- function(x)
  {
    if (!require(x,character.only = TRUE))
    {
      install.packages(x,dep=TRUE, repos='https://cloud.r-project.org/')
      if(!require(x,character.only = TRUE)) stop("Package not found")
    }
  }
pkgLoad("taxize")
message("")
message("Step 2: Import list of species")
taxon<-args[1]
mismatch<-read.table(paste0("MIDORI_",taxon,".ITIS_mismatch_sp.txt"),header=F,sep="\t", stringsAsFactors=FALSE)$V1
message("")
message("Step 3: Get EoL matches with classification()") # retrieve taxonomic hierarchy of the species in mismatch
ITISmismatch_class<-list()
for (i in mismatch) {
  SP_class<-classification(i, db="eol", return_id=FALSE, rows=1)
  ITISmismatch_class[i]<-SP_class
  rm(SP_class)
}
message("")
message("Step 4: Write successful EoL lookups to file")
library(rlist)
eol_matchedunlist<-ITISmismatch_class[unlist(lapply(ITISmismatch_class,function(x)is.character(x[[1]])))]
eol_matched <- list.clean(eol_matchedunlist, fun = is.character)
class(eol_matched) <- "classification"
eol_tab<-cbind(eol_matched)
write.table(cbind(rep("Eukaryota",nrow(eol_tab)),eol_tab$phylum,eol_tab$class,eol_tab$order,eol_tab$family,eol_tab$genus,eol_tab$species,rep("EoL_match",nrow(eol_tab)),eol_tab$query),paste0(taxon,"_EoL_matched_taxonomy.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
message("	Done. File is:")
paste0(taxon,"_EoL_matched_taxonomy.txt")
message("")

message("Step 5: Write failed lookups to file")
eol_unmatchedfalse <- list.clean(eol_matchedunlist, fun = is.data.frame)
eol_unmatchedNA <- ITISmismatch_class[unlist(lapply(ITISmismatch_class,function(x)is.logical(x[[1]])))]
eol_unmatched <- c(eol_unmatchedfalse,eol_unmatchedNA)
eol_missing <- names(eol_unmatched)
if((length(eol_missing))!=0){
  write.table(eol_missing,paste0(taxon,".EOL_missing_sp.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  message("	File is:")
  paste0(taxon,".EOL_missing_sp.txt")
} else {
  message("	No failed lookups to write out.")
}
message("")

message("End of get_taxonomy_mismatches_EOL.R")
q(save="no")
