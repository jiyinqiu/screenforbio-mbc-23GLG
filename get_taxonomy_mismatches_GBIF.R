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
mismatch<-read.table(paste0(taxon,".EOL_missing_sp.txt"),header=F,sep="\t", stringsAsFactors=FALSE)$V1
message("")
message("Step 3: Get GBIF matches with classification()") # retrieve taxonomic hierarchy of the species in mismatch
EOLmismatch_class<-list()
for (i in mismatch) {
  SP_class<-classification(i, db="gbif", return_id=FALSE, rows=1)
  EOLmismatch_class[i]<-SP_class
  rm(SP_class)
}
message("")
message("Step 4: Write successful gbif lookups to file")
library(rlist)
gbif_matchedunlist<-EOLmismatch_class[unlist(lapply(EOLmismatch_class,function(x)is.character(x[[1]])))]
gbif_matched <- list.clean(gbif_matchedunlist, fun = is.character)
class(gbif_matched) <- "classification"
gbif_tab<-cbind(gbif_matched)
write.table(cbind(rep("Eukaryota",nrow(gbif_tab)),gbif_tab$phylum,gbif_tab$class,gbif_tab$order,gbif_tab$family,gbif_tab$genus,gbif_tab$species,rep("GBIF_match",nrow(gbif_tab)),gbif_tab$query),paste0(taxon,"_GBIF_matched_taxonomy.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
message("	Done. File is:")
paste0(taxon,"_GBIF_matched_taxonomy.txt")
message("")

message("Step 5: Write failed lookups to file")
gbif_unmatchedfalse <- list.clean(gbif_matchedunlist, fun = is.data.frame)
gbif_unmatchedNA <- EOLmismatch_class[unlist(lapply(EOLmismatch_class,function(x)is.logical(x[[1]])))]
gbif_unmatched <- c(gbif_unmatchedfalse,gbif_unmatchedNA)
gbif_missing <- names(gbif_unmatched)
if((length(gbif_missing))!=0){
  write.table(gbif_missing,paste0(taxon,".GBIF_missing_sp.txt"),sep="\t",col.names=FALSE,row.names=FALSE,quote=FALSE)
  message("	File is:")
  paste0(taxon,".GBIF_missing_sp.txt")
} else {
  message("	No failed lookups to write out.")
}
message("")

message("End of get_taxonomy_mismatches_GBIF.R")
q(save="no")
