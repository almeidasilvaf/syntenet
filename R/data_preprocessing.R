
#' Check if input objects is ready for further analyses
#' 
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{annotation}.
#' @param annotation A GRangesList, CompressedGRangesList, or list of
#' GRanges with the annotation for the sequences in \strong{seq}. This list must
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
    
    valid <- valid_seq(seq) & valid_annot(annotation)

    check1 <- check_list_names(seq, annotation)
    check2 <- check_ngenes(seq, annotation)
    check3 <- check_gene_names(seq, annotation)
    
    return(TRUE)
}


#' Process sequence data
#'
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{annotation}.
#' @param annotation A GRangesList, CompressedGRangesList, or list of
#' GRanges with the annotation for the sequences in \strong{seq}. This list must
#' have names (not NULL), and names of each list element must match the names
#' of list elements in \strong{seq}.
#'
#' @return A list of 2 elements:
#' \describe{
#'   \item{seq}{The processed list of AAStringSet objects from \strong{seq}.}
#'   \item{annotation}{The processed GRangesList or CompressedGRangesList
#'   object from \strong{annotation}.}
#' }
#' 
#' @details 
#' This function processes the input sequences and annotation to:
#' \enumerate{
#'   \item Remove whitespace and anything after it in sequence names 
#'   (i.e., \code{names(seq[[x]])}, which is equivalent to FASTA headers), if
#'   there is any.
#'   \item Remove period followed by number at the end of sequence names, which
#'   typically indicates different isoforms of the same gene 
#'   (e.g., Arabidopsis thaliana's transcript AT1G01010.1, which belongs to
#'   gene AT1G01010).
#'   \item Add a unique species identifier to sequence names. The species 
#'   identifier consists of the first 3-5 strings of the element name.
#'   For instance, if the first element of the \strong{seq} list is named
#'   "Athaliana", each sequence in it will have an identifier "Atha_" added
#'   to the beginning of each gene name (e.g., Atha_AT1G01010).
#'   \item If sequences have an asterisk (*) representing stop codon, remove it.
#'   \item Add a unique species identifier (same as above) to 
#'   gene and chromosome names of each element of the \strong{annotation}
#'   GRangesList/CompressedGRangesList.
#'   \item Filter each element of the \strong{annotation} 
#'   GRangesList/CompressedGRangesList to keep only seqnames, ranges, and
#'   gene ID.
#' }
#' @importFrom Biostrings width subseq
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#' @rdname process_input
#' @examples 
#' data(annotation)
#' data(proteomes)
#' seq <- proteomes
#' clean_data <- process_input(seq, annotation)
process_input <- function(seq = NULL, annotation = NULL) {
    
    check <- check_input(seq, annotation)
    ab_len <- species_id_length(seq)
    s_abbrev <- substr(names(seq), 1, ab_len)
    a_abbrev <- substr(names(annotation), 1, ab_len)
    # Clean 'seq' object first
    clean_seqs <- lapply(seq_along(seq), function(x) {
        sequence <- seq[[x]]
        gene_names <- names(sequence)
        # 1. Remove everything after space
        gene_names <- gsub(" .*", "", gene_names)
        # 2. Remove period + number at the end of gene IDs
        gene_names <- gsub("\\.[0-9]$", "", gene_names)
        # 3. Add unique identifier in front of gene names
        gene_names <- paste(s_abbrev[x], gene_names, sep = "_")
        names(sequence) <- gene_names
        # 4. Remove "*" as stop codon (if any)
        last_aa <- sequence[[1]][Biostrings::width(sequence)[1]]
        if(identical(toString(last_aa), "*")) {
            sequence <- Biostrings::subseq(sequence, 1, width(sequence)-1)
        }
        return(sequence)
    })
    names(clean_seqs) <- names(seq)
    # Clean annotation object
    clean_annotation <- lapply(seq_along(annotation), function(x) {
        # 5. Add unique identifier in front of chromosome and gene names
        annot_df <- annotation[[x]][annotation[[x]]$type == "gene"]
        annot_df <- as.data.frame(annot_df)
        annot_df$seqnames <- paste(a_abbrev[x], annot_df$seqnames, sep = "_")
        if("gene_id" %in% colnames(annot_df)) {
            annot_df$gene <- paste(a_abbrev[x], annot_df$gene_id, sep = "_")
        } else {
            annot_df$gene <- paste(a_abbrev[x], annot_df$Name, sep = "_")
        }
        # 6. Keep only seqnames, ranges, and gene ID
        annot_df <- annot_df[, c("seqnames", "start", "end", "gene")]
        return(makeGRangesFromDataFrame(annot_df, keep.extra.columns = TRUE))
    })
    names(clean_annotation) <- names(annotation)
    flist <- list(seq = clean_seqs, annotation = clean_annotation)
    return(flist)
}


