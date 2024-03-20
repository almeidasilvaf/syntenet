
#' Binarize and transpose the phylogenomic profile matrix
#' 
#' @param profile_matrix A matrix with phylogenomic profiles obtained
#' with \code{phylogenomic_profile}.
#' 
#' @return A binary and transposed version of the profiles matrix.
#'
#' @export
#' @rdname binarize_and_transpose
#' @examples 
#' data(clusters)
#' profile_matrix <- phylogenomic_profile(clusters)
#' tmat <- binarize_and_transpose(profile_matrix)
binarize_and_transpose <- function(profile_matrix = NULL) {
    
    M <- profile_matrix

    # Make sure M is a numeric matrix of non-negative integer values
    if(!is.null(M)) {
        check_nm <- is.matrix(M) & is.numeric(M)
        check_nninteger <- as.vector(M >= 0 & M == round(M))
        bad_input <- any(c(check_nm, check_nninteger) == FALSE)
        if(bad_input) {
            stop("Input must be a matrix of non-negative integer values.")
        }
    }
    
    M[is.na(M)] <- 0
    M[M > 1] <- 1
    tM <- t(M)
    return(tM)
}


#' Save the transposed binary profiles matrix to a file in PHYLIP format
#'
#' @param transposed_profiles A binary and transposed profile matrix. 
#' The profile matrix can be obtained with \code{phylogenomic_profile()}.
#' @param outdir Path to output directory. By default, files are saved in
#' a temporary directory, so they will be deleted when the R session closes.
#' If you want to keep the files, specify a custom output directory.
#' 
#' @return Character specifying the path to the PHYLIP file.
#' @importFrom utils write.table
#' @export
#' @rdname profiles2phylip
#' @examples
#' data(clusters)
#' profile_matrix <- phylogenomic_profile(clusters)
#' tmat <- binarize_and_transpose(profile_matrix)
#' profiles2phylip(tmat)
profiles2phylip <- function(transposed_profiles = NULL, outdir = tempdir()) {
    
    dt <- format(Sys.time(), "%d_%b_%Y_%Hh%M")
    pat <- paste0("microsynteny_phylogeny_", dt)
    if(!dir.exists(outdir)) { 
        stop("Could not find directory specified in 'outdir'.") 
    }
    matrix_file <- file.path(outdir, paste0(pat, ".phy"))
    
    header <- paste(nrow(transposed_profiles), ncol(transposed_profiles))
    writeLines(header, matrix_file)
    write.table(transposed_profiles, file = matrix_file,
                col.names = FALSE, quote = FALSE, append = TRUE)
    
    return(matrix_file)
}

#' Infer microsynteny-based phylogeny with IQTREE
#'
#' @param transposed_profiles A binary and transposed profile matrix. 
#' The profile matrix can be obtained with \code{phylogenomic_profile()}.
#' @param bootr Numeric scalar with the number of bootstrap replicates.
#' Default: 1000.
#' @param alrtboot Numeric scalar with the number of replicates for
#' the SH-like approximate likelihood ratio test. Default: 1000.
#' @param threads Numeric scalar indicating the number of threads to use or
#' "AUTO", which allows IQTREE to automatically choose the best number 
#' of threads to use. Default: "AUTO".
#' @param model Substitution model to use. If you are unsure, pick the default.
#' Default: "MK+FO+R".
#' @param outdir Path to output directory. By default, files are saved in
#' a temporary directory, so they will be deleted when the R session closes.
#' If you want to keep the files, specify a custom output directory.
#' @param outgroup Name of outgroup clade to group the phylogeny. 
#' Default: NULL (unrooted phylogeny).
#' @param verbose Logical indicating if progress messages should be prompted.
#' Default: FALSE.
#' 
#' @return A character vector of paths to output files.
#' @export
#' @rdname infer_microsynteny_phylogeny
#' @examples 
#' data(clusters)
#' profile_matrix <- phylogenomic_profile(clusters)
#' tmat <- binarize_and_transpose(profile_matrix)
#' 
#' # Leave only some legumes and P. mume as an outgroup for testing purposes
#' included <- c("gma", "pvu", "vra", "van", "cca", "pmu")
#' tmat <- tmat[rownames(tmat) %in% included, ]
#' 
#' # Remove non-variable sites
#' tmat <- tmat[, colSums(tmat) != length(included)]
#' 
#' if(iqtree_is_installed()) {
#'     phylo <- infer_microsynteny_phylogeny(tmat, outgroup = "pmu", 
#'                                           threads = 1)
#' }
#'  
infer_microsynteny_phylogeny <- function(transposed_profiles = NULL,
                                         bootr = 1000, alrtboot = 1000,
                                         threads = "AUTO",
                                         model = "MK+FO+R",
                                         outdir = tempdir(),
                                         outgroup = NULL,
                                         verbose = FALSE) {
    
    if(!iqtree_is_installed()) { stop("Could not find IQ-TREE in PATH.") }
    
    # Create a file containing a header and the matrix below it
    matrix_file <- profiles2phylip(transposed_profiles, outdir)
    
    # Handle arguments
    stdout <- FALSE
    if(verbose) { stdout <- "" }
    
    root <- NULL
    if(!is.null(outgroup)) {
        root <- paste0("-o ", outgroup)
    }
    
    # Run IQ-TREE
    if(iqtree_version() == 1) { # IQ-TREE v1
        iqtree_args <- c(
            "-s ", matrix_file, "-bb", bootr, "-alrt", alrtboot, "-nt", threads, root, "-m", model, "-st MORPH -redo"
        )
        iqtree <- system2("iqtree", args = iqtree_args, stdout = stdout)
    } else { # IQ-TREE v2
        iqtree_args <- c(
            "-s ", matrix_file, 
            "-B", bootr, 
            "-alrt", alrtboot, 
            "-T", threads, 
            root,
            "-m", model,
            "-st MORPH -redo"
        )
        iqtree <- system2("iqtree2", args = iqtree_args, stdout = stdout)
    }
    
    paths <- list.files(
        path = outdir, pattern = basename(matrix_file), full.names = TRUE
    )
    return(paths)
}
