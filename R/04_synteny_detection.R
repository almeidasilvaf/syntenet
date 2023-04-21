
#' Create a data frame of pairwise comparisons to make in \code{run_diamond}
#' 
#' @param species Character vector of species names.
#' @param compare Character scalar indicating which comparisons
#' should be made when running DIAMOND. 
#' Possible modes are "all" (all-vs-all comparisons), 
#' "intraspecies" (intraspecies comparisons only), or
#' "interspecies" (interspecies comparisons only). Alternatively, users can
#' pass a 2-column data frame as input with the names of species to be 
#' compared.
#'
#' @return A 2-column data frame with columns \strong{sp1} and \strong{sp2}
#' representing the pairwise species comparisons to be made.
#' @importFrom methods is
#' @noRd
#' @examples 
#' species <- c("Olucimarinus", "OspRCC809")
#' get_comp(species)
get_comp <- function(species, compare = "all") {

    # Create a data frame of all possible pairwise comparisons
    comb_df <- expand.grid(species, species, stringsAsFactors = FALSE)
    names(comb_df) <- c("sp1", "sp2")
    
    if(is(compare, "data.frame")) {
        if(ncol(compare) != 2) {
            stop("The data frame of species comparisons must have 2 columns.")
        }
        comb_df <- compare
        names(comb_df) <- c("sp1", "sp2")
        
        # Check if any species in input data frame is not in 'seq'
        sp <- unique(c(comb_df$sp1, comb_df$sp2))
        name_is_missing <- any(!sp %in% species)
        if(name_is_missing) {
            stop("All species names in 'compare' must be in 'seq'.")
        }
        
    } else if(compare == "intraspecies") {
        comb_df <- comb_df[comb_df$sp1 == comb_df$sp2, ]
        
    } else if(compare == "interspecies") {
        comb_df <- comb_df[comb_df$sp1 != comb_df$sp2, ]
        
        
    } else if(compare == "all") {
        comb_df <- comb_df
        
    } else {
        stop("Invalid argument passed to parameter 'compare'.")
    }
    return(comb_df)
}

#' Get edges (anchor pairs) from .collinearity files obtained with MCScan
#' 
#' @param collinearity_paths Character vector of paths to .collinearity files.
#' 
#' @return A data frame with variables:
#' \describe{
#'   \item{Block}{Numeric, synteny block ID}
#'   \item{Anchor1}{Character, gene ID of anchor 1.}
#'   \item{Anchor2}{Character, gene ID of anchor 2.}
#' }
#' 
#' @importFrom utils read.table
#' @noRd
edges_from_collinearity <- function(collinearity_paths = NULL) {
    
    edges <- lapply(seq_along(collinearity_paths), function(x) {
        lines <- readLines(collinearity_paths[x])
        nlines <- length(lines[!startsWith(lines, "#")])
        
        df <- NULL
        if(nlines > 0) {
            df <- read.table(collinearity_paths[x], sep = "\t", comment.char = "#")
            df <- df[, c(1, 2, 3)]
            names(df) <- c("Block", "Anchor1", "Anchor2")
            df$Block <- as.numeric(gsub("-.*", "", df$Block))
            
        }
        return(df)
    })
    
    edges <- Reduce(rbind, edges)
    return(edges)
}

#' Get synteny blocks from .collinearity files obtained with MCScan
#' 
#' @param collinearity_paths Character vector of paths to .collinearity files.
#' 
#' @return A data frame with variables:
#' \describe{
#'   \item{Block}{Numeric, synteny block ID}
#'   \item{Block_score}{Numeric, score of synteny block.}
#'   \item{Chr}{Character, query and target chromosome of the synteny
#'              block formatted as "<querychr>&<targetchr>".}
#'   \item{Orientation}{Character, the orientation of genes within blocks,
#'                      with "plus" indicating that genes are in the same 
#'                      direction, and "minus" indicating that genes are 
#'                      in opposite directions.}
#' }
#' 
#' @noRd
blocks_from_collinearity <- function(collinearity_paths = NULL) {
    
    blocks <- lapply(seq_along(collinearity_paths), function(x) {
        lines <- readLines(collinearity_paths[x])
        nlines <- length(lines[!startsWith(lines, "#")])
        
        block_df <- NULL
        if(nlines > 0) {
            
            # Extract synteny block information
            brows <- lines[grep("^## Alignment", lines)]
            binfo <- strsplit(gsub("## Alignment ", "", brows), split = " ")
            
            # Create a data frame with synteny block information
            block_df <- as.data.frame(do.call(rbind, binfo))[, c(1, 2, 5, 6)]
            names(block_df) <- c("Block", "Block_score", "Chr", "Orientation")
            block_df$Block <- as.numeric(gsub(":.*", "", block_df$Block))
            block_df$Block_score <- as.numeric(
                gsub("score=", "", block_df$Block_score)
            )
        }
        return(block_df)
    })
    
    blocks <- Reduce(rbind, blocks)
    return(blocks)
}


#' Parse .collinearity files obtained with MCScan
#' 
#' The .collinearity files can be obtained with \code{intraspecies_synteny}
#' and \code{interspecies_synteny}, which execute a native version of the
#' MCScan algorithm.
#'
#' @param collinearity_paths Character vector of paths to .collinearity files.
#' @param as Character specifying what to extract. 
#' One of "anchors" (default), "blocks", or "all".
#'
#' @return If \strong{as} is "anchors", a data frame with variables "Anchor1",
#' and "Anchor2". If \strong{as} is "blocks", a data frame with 
#' variables "Block", "Block_score", "Chr", and "Orientation". 
#' If \strong{as} is "all", a data frame with all aforementioned variables,
#' which indicate:
#' \describe{
#'   \item{Block}{Numeric, synteny block ID}
#'   \item{Block_score}{Numeric, score of synteny block.}
#'   \item{Chr}{Character, query and target chromosome of the synteny
#'              block formatted as "<querychr>&<targetchr>".}
#'   \item{Orientation}{Character, the orientation of genes within blocks,
#'                      with "plus" indicating that genes are in the same 
#'                      direction, and "minus" indicating that genes are 
#'                      in opposite directions.}
#'   \item{Anchor1}{Character, gene ID of anchor 1.}
#'   \item{Anchor2}{Character, gene ID of anchor 2.}
#' }
#'
#' @importFrom utils read.table
#' @export
#' @rdname parse_collinearity
#' @examples 
#' collinearity_paths <- system.file(
#'     "extdata", "Scerevisiae.collinearity", package = "syntenet"
#' )
#' net <- parse_collinearity(collinearity_paths)
parse_collinearity <- function(
        collinearity_paths = NULL, as = "anchors"
) {
    
    # Add clean names to file paths (basename)
    fname <- gsub("\\.collinearity", "", basename(collinearity_paths))
    names(collinearity_paths) <- fname
    
    # Get parsed files
    if(as == "anchors") {
        parsed <- edges_from_collinearity(collinearity_paths)
        parsed <- parsed[, c("Anchor1", "Anchor2")]
        
    } else if(as == "blocks") {
        parsed <- blocks_from_collinearity(collinearity_paths)
        
    } else if(as == "all") {
        parsed <- merge(
            blocks_from_collinearity(collinearity_paths),
            edges_from_collinearity(collinearity_paths)
        )
    } else {
        
        stop("Invalid argument passed to 'as' parameter.")
    }
    
    return(parsed)
}


#' Infer synteny network
#' 
#' @param blast_list A list of data frames, each data frame having 
#' the tabular output of BLASTp or similar programs, such as DIAMOND. 
#' This is the output of the function \code{run_diamond()}. If you performed
#' pairwise comparisons on the command line, you can read the tabular output
#' as data frames and combine them in a list. List names must be have species
#' names separated by underscore. For instance, if the first list element is
#' a data frame containing the comparison of speciesA (query) 
#' against speciesB (database), its name must be "speciesA_speciesB".
#' @param annotation A processed GRangesList, CompressedGRangesList, or
#' list of GRanges as returned by \code{process_input()}. 
#' @param outdir Path to the output directory. Default: tempdir().
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param is_pairwise specify if only pairwise blocks should be reported
#' Default: TRUE
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param ... Any additional arguments to
#' `mcscanx`.
#' 
#' @return A network represented as an edge list.
#'
#' @export
#' @rdname infer_syntenet
#' @importFrom utils write.table
#' @examples
#' # Load data
#' data(proteomes)
#' data(annotation)
#' data(blast_list)
#'
#' # Create processed annotation list
#' annotation <- process_input(proteomes, annotation)$annotation
#'
#' # Infer the synteny network
#' net <- infer_syntenet(blast_list, annotation)
infer_syntenet <- function(
        blast_list = NULL, 
        annotation = NULL, 
        outdir = tempdir(),
        anchors = 5, max_gaps = 25, is_pairwise = TRUE, 
        verbose = FALSE, ...
) {

    valid <- valid_blast(blast_list) & valid_annot(annotation)
    
    # Separate intra from interspecies
    species <- names(annotation)
    equal_comp <- vapply(species, function(x) paste0(x, "_", x), character(1))
    idx_equal <- which(names(blast_list) %in% equal_comp)
    
    #---- 1) Intraspecies synteny detection------------------------------------
    intra_dir <- file.path(outdir, "intraspecies_synteny")
    if(!dir.exists(intra_dir)) { dir.create(intra_dir, recursive = TRUE) }
    blast_intra <- blast_list[idx_equal]
    
    intraspecies <- intraspecies_synteny(
        blast_intra = blast_intra, 
        annotation = annotation,
        intra_dir = intra_dir, 
        anchors = anchors, 
        max_gaps = max_gaps, 
        is_pairwise = is_pairwise,
        verbose = verbose, 
        ...
    )
    
    #---- 2) Interspecies synteny detection------------------------------------
    inter_dir <- file.path(outdir, "interspecies_synteny")
    if(!dir.exists(inter_dir)) { dir.create(inter_dir, recursive = TRUE) }
    blast_inter <- blast_list[-idx_equal]
    
    interspecies <- interspecies_synteny(
        blast_inter = blast_inter, 
        annotation = annotation, 
        inter_dir = inter_dir,
        anchors = anchors, 
        max_gaps = max_gaps, 
        is_pairwise = is_pairwise,
        verbose = verbose, 
        ...
    )
    
    # Create edge list
    syn_files <- c(intraspecies, interspecies)
    edges <- parse_collinearity(syn_files)
    return(edges)
}
