
#' Read GFF/GTF files in a directory as a GRangesList object
#'
#' @param gff_dir Character indicating the path to the directory containing
#' GFF/GTF files.
#'
#' @return A GRangesList object, where each element represents a different
#' GFF/GTF file.
#' 
#' @export
#' @rdname gff2GRangesList
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRangesList
#' @examples 
#' gff_dir <- system.file("extdata", "annotation", package = "syntenet")
#' grangeslist <- gff2GRangesList(gff_dir)
gff2GRangesList <- function(gff_dir) {

    # Get full path to GFF/GFF3/GTF files    
    files <- list.files(gff_dir, pattern = "\\.gff|\\.gtf", full.names = TRUE)

    # Import files
    grangeslist <- lapply(files, rtracklayer::import)
    grangeslist <- GenomicRanges::GRangesList(grangeslist)
    
    # Name elements based on file names
    names(grangeslist) <- gsub("\\.gff.*|\\.gtf.*", "", basename(files))
    
    return(grangeslist)
}


#' Read FASTA files in a directory as a list of AAStringSet objects
#' 
#' @param fasta_dir Character indicating the path to the directory containing
#' FASTA files.
#' 
#' @return A list of AAStringSet objects, where each element represents a 
#' different FASTA file.
#' 
#' @export
#' @rdname fasta2AAStringSetlist
#' @importFrom Biostrings readAAStringSet 
#' @examples
#' fasta_dir <- system.file("extdata", "sequences", package = "syntenet")
#' aastringsetlist <- fasta2AAStringSetlist(fasta_dir)
fasta2AAStringSetlist <- function(fasta_dir) {
    
    # Get full path to FASTA files    
    files <- list.files(
        fasta_dir, pattern = "\\.fa|\\.fna", full.names = TRUE
    )
    
    # Read FASTA files as a list of `AAStringSet` objects
    seqs <- lapply(files, Biostrings::readAAStringSet)
    
    # Name list elements
    names(seqs) <- gsub("\\.fa.*|\\.fna.*", "", basename(files))
    
    return(seqs)
}
