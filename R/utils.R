
#' Check if the names of list of sequences and annotations match.
#' 
#'
#' @param seq A list of AAStringSet objects.
#' @param annotation A GRangesList or CompressedGRangesList object with
#' the annotation for the sequences in \strong{seq}.
#' 
#' @return TRUE if the objects pass the check.
#' @noRd
#' @examples 
#' data(annotation)
#' data(proteomes)
#' seq <- proteomes
#' check_list_names(seq, annotation)
check_list_names <- function(seq = NULL, annotation = NULL) {
    
    annot_names <- names(annotation)
    seq_names <- names(seq)
    n_match <- annot_names %in% seq_names
    
    if(is.null(annot_names) | is.null(seq_names)) {
        stop("List-like arguments 'seq' and 'annotation' must have names.")
    } else if(any(isFALSE(n_match))) {
        stop("Names of list elements in 'seq' and 'annotation' must match.")
    } else {
        check <- TRUE
    }
    return(check)
}

#' Check if the number of sequences is less than the number of genes
#'
#' @param seq A list of AAStringSet objects.
#' @param annotation A GRangesList or CompressedGRangesList object 
#' with the annotation for the sequences in \strong{seq}.
#'
#' @return TRUE if the objects pass the check.
#' @noRd 
#' @examples 
#' data(proteomes)
#' data(annotation)
#' seq <- proteomes
#' check_ngenes(seq, annotation)
check_ngenes <- function(seq = NULL, annotation = NULL) {
    
    # Data frame of species and gene count based on annotation
    gene_count <- Reduce(rbind, lapply(seq_along(annotation), function(x) {
        count <- length(annotation[[x]][annotation[[x]]$type == "gene"])
        count_df <- data.frame(
            Species = names(annotation)[x],
            Genes = count
        )
        return(count_df)
    }))
    
    # Data frame of species and gene count based on sequences
    seq_count <- Reduce(rbind, lapply(seq_along(seq), function(x) {
        count <- length(seq[[x]])
        count_df <- data.frame(
            Species = names(annotation)[x],
            Seqs = count
        )
        return(count_df)
    }))
    
    # Check if number of sequences is <= gene count (accounting for ncRNAs)
    counts <- merge(gene_count, seq_count, by = "Species")
    check_count <- counts$Seqs <= counts$Genes
    idx_error <- which(check_count == FALSE)
    if(length(idx_error) != 0) {
        name <- counts$Species[idx_error]
        name <- paste0(seq_along(name), ". ", name)
        name <- paste0(name, collapse = "\n")
        stop("Number of sequences in greater than the number of genes for:\n",
             name)
    } 
    return(TRUE)
}


#' Check if sequence names match gene names from annotation
#'
#' @param seq A list of AAStringSet objects.
#' @param annotation A GRangesList or CompressedGRangesList object 
#' with the annotation for the sequences in \strong{seq}.
#'
#' @return TRUE if the objects pass the check.
#' @noRd 
#' @importFrom GenomicRanges mcols
#' @examples
#' data(annotation)
#' data(proteomes)
#' seq <- proteomes
#' check_gene_names(seq, annotation)
check_gene_names <- function(seq = NULL, annotation = NULL) {
    
    seq_names <- lapply(seq, names)
    gene_names <- lapply(annotation, function(x) {
        cols <- names(GenomicRanges::mcols(x))
        if("gene_id" %in% cols) {
            nam <- unique(x$gene_id)
        } else if("Name" %in% cols) {
            nam <- unique(x$Name)
        } else {
            stop("Could not find 'gene_id' or 'Name' column in GRanges.")
        }
        return(nam)
    })

    # Check if names in `seq` match gene names in `annotation`
    check_names <- lapply(seq_along(seq_names), function(x) {
        c <- seq_names[[x]] %in% gene_names[[x]]
        c <- any(c == FALSE)
        return(c)
    })
    
    idx_error <- which(check_names == TRUE) # TRUE means error
    if(length(idx_error) != 0) {
        name <- names(seq_names)[idx_error]
        name <- paste0(seq_along(name), ". ", name)
        name <- paste0(name, collapse = "\n")
        stop("Sequence names in 'seq' do not match gene names in 'annotation' for:\n",
             name)
    }
    return(TRUE)
}


#' Pick best length for unique species identifiers
#'
#' @param input_list A list of AAStringSet objects or a 
#' GRangesList/CompressedGRangesList object.
#' 
#' @return Numeric scalar with the length of the 
#' species ID (either 3, 4, or 5).
#' @noRd
species_id_length <- function(input_list = NULL) {
    
    nam <- names(input_list)
    nam <- gsub(" ", "_", nam)
    abbrev <- unlist(lapply(nam, substr, 1, 3)) # try with 3
    count <- table(abbrev)
    n <- 3
    if(any(count > 1)) {
        abbrev <- unlist(lapply(nam, substr, 1, 4)) # try with 4
        count <- table(abbrev)
        n <- 4
        if(any(count > 1)) {
            abbrev <- unlist(lapply(nam, substr, 1, 5)) # try with 5
            count <- table(abbrev)
            n <- 5
        }
    }
    return(n)
}


#' Wrapper to check if command is found in PATH
#'
#' @param cmd Command to test.
#' @param args Arguments for command.
#'
#' @author Fabricio Almeida-Silva
#' @return Logical indicating whether the command is in PATH or not.
#' @noRd
is_valid <- function(cmd = NULL, args = NULL) {
    found <- tryCatch(
        system2(cmd, args = args, stdout = FALSE, stderr = FALSE),
        error = function(e) return(FALSE),
        warning = function(w) return(FALSE)
    )
    if(!isFALSE(found)) {
        found <- TRUE
    }
    return(found)
}


#' Check if DIAMOND is installed
#'
#' @return Logical indicating whether DIAMOND is installed or not.
#' @export
#' @rdname diamond_is_installed
#' @examples
#' diamond_is_installed()
diamond_is_installed <- function() {
    valid <- is_valid(cmd = "diamond", args = "-h")
    return(valid)
}
