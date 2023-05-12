

#' Collapse protein IDs into gene IDs in sequence names of AAStringSet objects
#' 
#' This function can be used if the sequence names of the AAStringSet objects
#' contain protein IDs instead of gene IDs (what syntenet requires)
#'
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{protein2gene}.
#' @param protein2gene A list of 2-column data frames containing
#' protein-to-gene ID correspondences, where the first column contains 
#' protein IDs, and the second column contains gene IDs. Names of list elements
#' must match names of \strong{seq}.
#' 
#' @return A list of AAStringSet objects as in \strong{seq}, but with
#' protein IDs replaced with gene IDs.
#' 
#' @details
#' For each species, this function will replace the protein IDs in sequence
#' names with gene IDs using the protein-to-gene correspondence table in
#' \strong{protein2gene}. After replacing protein IDs with gene IDs, if
#' there are multiple sequences with the same gene ID (indicating different
#' isoforms of the same gene), only the longest sequence is kept, so that
#' the number of sequences is not greater than the number of genes.
#' 
#' @importFrom Biostrings width
#' @export
#' @rdname collapse_protein_ids
#' @examples
#' # Load data
#' seq_path <- system.file(
#'     "extdata", "RefSeq_parsing_example", package = "syntenet"
#' )
#' seq <- fasta2AAStringSetlist(seq_path)
#' annot <- gff2GRangesList(seq_path)
#' 
#' # Clean sequence names
#' names(seq$Aalosa) <- gsub(" .*", "", names(seq$Aalosa))
#' 
#' # Create a correspondence data frame
#' cor_df <- as.data.frame(annot$Aalosa[annot$Aalosa$type == "CDS", ])
#' cor_df <- cor_df[, c("Name", "gene")]
#' 
#' # Create a list of correspondence data frames
#' protein2gene <- list(Aalosa = cor_df)
#' 
#' # Collapse IDs
#' new_seqs <- collapse_protein_ids(seq, protein2gene)
collapse_protein_ids <- function(seq, protein2gene = NULL) {
    
    if(!is(protein2gene, "list") | !is(protein2gene[[1]], "data.frame")) {
        stop("`protein2gene` must be a list of data frames.")
    }
    
    if(!all(names(seq) %in% names(protein2gene))) {
        stop("All species in `seq` must be in `protein2gene`.")
    }
    
    new_seqs <- lapply(seq_along(seq), function(x) {
        
        species <- names(seq)[x]
        
        nseq <- seq[[species]]
        pgene <- protein2gene[[species]]
        
        ## Get indices of matches
        idx <- unique(match(
            names(nseq), 
            pgene[, 1]
        ))
        if(length(idx) <= 1) {
            
            m <- paste0(
                "None of the protein IDs in column 1 of `protein2gene[['",
                species, "']]` match sequence names in `seq[['", 
                species, "']]`."
                
            )
            stop(m)
        }
        
        ## Replace names
        names(nseq) <- pgene[, 2][idx]
        
        # Keep only longest protein for each gene
        nseq <- nseq[order(Biostrings::width(nseq), decreasing = TRUE), ]
        nseq <- nseq[!duplicated(names(nseq)), ]
        
        return(nseq)
    })
    names(new_seqs) <- names(seq)
    
    return(new_seqs)
    
}



#' Check if input objects are ready for further analyses
#' 
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{annotation}.
#' @param annotation A GRangesList, CompressedGRangesList, or list of
#' GRanges with the annotation for the sequences in \strong{seq}. This list must
#' have names (not NULL), and names of each list element must match the names
#' of list elements in \strong{seq}.
#' @param gene_field Character, name of the column in the GRanges objects
#' that contains gene IDs. Default: "gene_id".
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
check_input <- function(seq = NULL, annotation = NULL, gene_field = "gene_id") {
    
    valid <- valid_seq(seq) & valid_annot(annotation)

    check1 <- check_list_names(seq, annotation)
    check2 <- check_ngenes(seq, annotation)
    check3 <- check_gene_names(seq, annotation, gene_field = gene_field)
    
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
#' @param gene_field Character, name of the column in the GRanges objects
#' that contains gene IDs. Default: "gene_id".
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
process_input <- function(seq = NULL, annotation = NULL, 
                          gene_field = "gene_id") {
    
    check <- check_input(seq, annotation, gene_field = gene_field)
    
    seq_ids <- create_species_id_table(names(seq))$species_id
    annot_ids <- create_species_id_table(names(annotation))$species_id
    
    # Clean 'seq' object first
    clean_seqs <- lapply(seq_along(seq), function(x) {
        sequence <- seq[[x]]
        gene_names <- names(sequence)
        
        ## 1. Remove everything after space
        gene_names <- gsub(" .*", "", gene_names)
        
        ## 2. Remove period + number at the end of gene IDs
        gene_names <- gsub("\\.[0-9]$", "", gene_names)
        
        ## 3. Add unique identifier in front of gene names
        gene_names <- paste(seq_ids[x], gene_names, sep = "_")
        names(sequence) <- gene_names
        
        ## 4. Remove "*" as stop codon (if any)
        last_aa <- sequence[[1]][Biostrings::width(sequence)[1]]
        if(identical(toString(last_aa), "*")) {
            sequence <- Biostrings::subseq(sequence, 1, width(sequence)-1)
        }
        return(sequence)
    })
    names(clean_seqs) <- names(seq)
    
    # Clean annotation object
    clean_annotation <- lapply(seq_along(annotation), function(x) {
        
        ## 5. Add unique identifier in front of chromosome and gene names
        annot_df <- annotation[[x]][annotation[[x]]$type == "gene"]
        annot_df <- as.data.frame(annot_df)
        annot_df$seqnames <- paste(annot_ids[x], annot_df$seqnames, sep = "_")
        annot_df$gene <- paste(annot_ids[x], annot_df[[gene_field]], sep = "_")

        # 6. Keep only seqnames, ranges, and gene ID
        annot_df <- annot_df[, c("seqnames", "start", "end", "strand", "gene")]
        return(makeGRangesFromDataFrame(annot_df, keep.extra.columns = TRUE))
    })
    names(clean_annotation) <- names(annotation)
    flist <- list(seq = clean_seqs, annotation = clean_annotation)
    return(flist)
}


