
#' Check if the names of list of sequences and annotations match.
#' 
#'
#' @param seq A list of AAStringSet objects.
#' @param annotation A GRangesList, CompressedGRangesList, or list of
#' GRanges with the annotation for the sequences in \strong{seq}.
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
    } else if(any(n_match == FALSE)) {
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
    
    for(n in seq(3, 5)) { 
        abbrev <- lapply(nam, substr, 1, n)
        count <- table(unlist(abbrev))
        if (!any(count > 1)) break
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
    valid <- is_valid(cmd = "diamond", args = "help")
    return(valid)
}


#' Check if IQTREE is installed
#'
#' @return Logical indicating whether IQTREE is installed or not.
#' @export
#' @rdname iqtree_is_installed
#' @examples
#' iqtree_is_installed()
iqtree_is_installed <- function() {
    valid <- is_valid(cmd = "iqtree", args = "-h")
    return(valid)
}


#' Generate custom color palette
#'
#' @param pal Numeric specifying palette number, from 1 to 3.
#'
#' @return Character vector of custom color palette with 20 colors
#' @noRd
custom_palette <- function(pal = 1) {
    pal1 <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
              "#9467BDFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF",
              "#BCBD22FF", "#17BECFFF", "#AEC7E8FF", "#FFBB78FF",
              "#98DF8AFF", "#FF9896FF", "#C5B0D5FF", "#C49C94FF",
              "#F7B6D2FF", "#C7C7C7FF", "#DBDB8DFF", "#9EDAE5FF")
    
    pal2 <- c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF",
              "#636363FF", "#6BAED6FF", "#FD8D3CFF", "#74C476FF",
              "#9E9AC8FF", "#969696FF", "#9ECAE1FF", "#FDAE6BFF",
              "#A1D99BFF", "#BCBDDCFF", "#BDBDBDFF", "#C6DBEFFF",
              "#FDD0A2FF", "#C7E9C0FF", "#DADAEBFF", "#D9D9D9FF")
    
    pal3 <- c("#393B79FF", "#637939FF", "#8C6D31FF", "#843C39FF",
              "#7B4173FF", "#5254A3FF", "#8CA252FF", "#BD9E39FF",
              "#AD494AFF", "#A55194FF", "#6B6ECFFF", "#B5CF6BFF",
              "#E7BA52FF", "#D6616BFF", "#CE6DBDFF", "#9C9EDEFF",
              "#CEDB9CFF", "#E7CB94FF", "#E7969CFF", "#DE9ED6FF")
    
    l <- list(pal1, pal2, pal3)
    l_final <- l[[pal]]
    return(l_final)
}


#' Wrapper to handle gene color annotation in heatmap
#'
#' @param species_annotation A 2-column data frame with species IDs in 
#' the first column (same as column names of profile matrix), and species
#' annotation (e.g., higher-level taxonomic information) in the second column.
#'
#' @return A list containing a data frame of species annotation and a list with
#' a named character vector of colors for each annotation category.
#' @noRd
#' @examples 
#' species_order <- c(
#'     "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
#'     "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
#'     "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
#' ) 
#' species_annotation <- data.frame(
#'    Species = species_order,
#'    Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
#'               "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3), 
#'               "Linaceae", "Salicaceae")
#')
#' species_colors(species_annotation)
species_colors <- function(species_annotation) {
    
    # First column to rownames
    rownames(species_annotation) <- species_annotation[,1]
    species_annotation[,1] <- NULL
    
    a_col <- colnames(species_annotation)[1] # Name of annotation column
    n_categories <- length(unique(species_annotation[[a_col]]))
    
    # List of named character vector with annotation colors for each category
    annotation_colors <- list(
        custom_palette(2)[seq_len(n_categories)]
    )
    names(annotation_colors) <- a_col
    names(annotation_colors[[a_col]]) <- unique(species_annotation[[a_col]])
    
    # Save results in a list
    results <- list(
        species_annotation = species_annotation,
        annotation_colors = annotation_colors
    )
    return(results)
}


#' Wrapper to check if input sequences have the expected class
#'
#' @param seq A list of AAStringSet objects, each list element containing
#' protein sequences for a given species. This list must have names 
#' (not NULL), and names of each list element must match the names of
#' list elements in \strong{annotation}.
#' 
#' @return TRUE if the input is valid.
#' @noRd
#' @importFrom methods is
valid_seq <- function(seq) {
    is_valid <- FALSE
    
    if(is(seq, "list") & is(seq[[1]], "AAStringSet")) {
        is_valid <- TRUE
    }
    
    if(!is_valid) { 
        stop("Input sequences must be in a list of AAStringSet objects.") 
    }
    
    return(is_valid)
}


#' Wrapper to check if input annotation has the expected class
#'
#' @param annotation A GRangesList, CompressedGRangesList, or list of
#' GRanges with the annotation for the sequences in \strong{seq}. This list must
#' have names (not NULL), and names of each list element must match the names
#' of list elements in \strong{seq}.
#' 
#' @return TRUE if the input is valid.
#' @noRd
#' @importFrom methods is
valid_annot <- function(annot) {

    if(is(annot, "list") | is(annot, "GRangesList") | 
       is(annot, "CompressedGRangesList")) {
        is_valid <- TRUE
    } else {
        stop("Input annotation must be in a GRangesList, CompressedGRangesList,
             or list of GRanges.") 
    }
    
    return(is_valid)
}


#' Wrapper to check if input BLAST list has the expected class
#'
#' @param blast_list A list of data frames, each data frame having 
#' the tabular output of BLASTp or similar programs, such as DIAMOND.
#' 
#' @return TRUE if the input is valid.
#' @noRd
#' @importFrom methods is
valid_blast <- function(blast_list) {

    if(is(blast_list, "list") & is(blast_list[[1]], "data.frame")) {
        is_valid <- TRUE
    } else {
        stop("BLAST list must be a list of data frames.")
    }
    
    return(is_valid)
}

