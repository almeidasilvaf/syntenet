
#' Check if input objects is ready for further analyses
#' 
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{annotation}.
#' @param annotation A GRangesList or CompressedGRangesList object 
#' with the annotation for the sequences in \strong{seq}. This list must
#' have names (not NULL), and names of each list element must match the names
#' of list elements in \strong{seq}.
#'
#' @details 
#' This function checks the input data for 3 required conditions:
#' \enumerate{
#'   \item Names of \strong{seq} list (i.e., \code{names(seq)}) match 
#'   the names of \strong{annotation} GRangesList/CompressedGRangesList
#'   (i.e., \code{names(annotation)})
#'   \item For each species (list elements), the number of sequences
#'   in \strong{seq} is not greater than the number of genes 
#'   in \strong{annotation}. This is a way to ensure users do not input
#'   the translated sequences for multiple isoforms of the same gene (generated
#'   by alternative splicing). Ideally, the number of sequences in \strong{seq} 
#'   should be equal to the number of genes in \strong{annotation}, but
#'   this may not always stand true because of non-protein-coding genes.
#'   \item For each species, sequence names (i.e., \code{names(seq[[x]])},
#'   equivalent to FASTA headers) match gene names in \strong{annotation}.
#' }
#' @return TRUE if the objects pass the check.
#' @export
#' @rdname check_input 
#' @examples 
#' data(annotation) 
#' data(proteomes)
#' check_input(proteomes, annotation)
check_input <- function(seq = NULL, annotation = NULL) {
    
    check1 <- check_list_names()
    check2 <- check_ngenes()
    check3 <- check_gene_names()
    
    return(TRUE)
}


#' Process sequence data
#'
#' @param seq A list of AAStringSet objects.
#'
#' @examples 
#' data(proteomes)
process_input <- function(seq = NULL) {
    
}