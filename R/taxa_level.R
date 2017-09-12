#'Taxa level collating
#'
#'This function takes a physeq object and returns a physeq object with taxa at a specified taxonomic level.
#'It extends the \link[phyloseq]{tax_glom}  to include names of corresponding taxa 
#'at a specified taxonomy level and  updating tree tip labels accordingly.
#'
#' @param physeq (Required). A \link[phyloseq]{phyloseq} object containing merged information of abundance,
#'        taxonomic assignment, sample data including the measured variables and categorical information
#'        of the samples, and / or phylogenetic tree if available.
#' @param which_level (Required). Character string specifying taxonomic level.
#' @return physeq object at specified taxonomic level.
#' 
#' @import phyloseq
#' @import phangorn
#' 
#' @examples 
#' data(pitlatrine)
#' physeq<-data(physeq)
#' physeq <- taxa_level(physeq = physeq,which_level = "Family")
#'  @export taxa_level
#' 

taxa_level <- function(physeq,which_level,include.tree=F){
  
  tax_trim <- tax_glom(physeq,which_level)
  
  abund_table <- otu_table(tax_trim)
  meta_table <- sample_data(tax_trim)
  tax <- tax_table(tax_trim)
  
  new_names<-unique(tax[,which_level]) #get taxa names at specified level and make them new rownames
  
  colnames(abund_table) <- new_names #colnames to match with taxa table
  
  if(!include.tree){
    physeq<-merge_phyloseq(phyloseq(abund_table), meta_table)
  }
  else if(include.tree){
    rownames(tax) <- new_names
    ptree <- phy_tree(tax_trim)
    ptree$tip.label <- new_names #change tip labels of tree to match taxa and otu table
    physeq<-merge_phyloseq(abund_table,tax,meta_table,ptree)
  }
  return(physeq)
}
