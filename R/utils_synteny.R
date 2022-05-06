
#' Detect intraspecies synteny
#'
#' @param blast_intra A list of BLAST tables for intraspecies comparisons.
#' @param intra_dir Output directory.
#' @param annot_dfs A list of annotation data frames.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param ... Any additional arguments to
#' `mcscanx`.
#'
#' @return Paths to .collinearity files.
#' @noRd
intraspecies_synteny <- function(blast_intra = NULL, intra_dir = NULL,
                                 annot_dfs = NULL, anchors = 5, 
                                 max_gaps = 25, verbose = FALSE, ...) {
    
    intrasyn <- unlist(lapply(blast_intra, function(x) {
        sp <- gsub("_.*", "", x[1, 1])
        
        # Write files
        blast_file <- file.path(intra_dir, paste0(sp, ".blast"))
        write.table(x, file = blast_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        gff_file <- file.path(intra_dir, paste0(sp, ".gff"))
        write.table(annot_dfs[[sp]], file = gff_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        # Detect synteny
        input <- file.path(intra_dir, sp)
        if(verbose) { verbose <- "" }
        #syn_args <- c(input, "-s ", anchors, "-m ", max_gaps)
        #syn <- system2("MCScanX", args = syn_args, stdout = verbose)
        rcpp_mcscanx_file(blast_file=blast_file, gff_file=gff_file, prefix=sp,
            outdir=intra_dir, match_size=anchors, max_gaps=max_gaps,
            is_pairwise=TRUE, in_synteny=1, ...)
        
        # Delete intermediate files
        unlink(c(blast_file, gff_file))
        
        synf <- file.path(intra_dir, paste0(sp, ".collinearity"))
        
        return(synf)
    }))
    intrasyn <- intrasyn[!is.null(intrasyn)]
    return(intrasyn)
}


#' Detect interspecies synteny
#'
#' @param blast_inter A list of BLAST tables for interspecies comparisons.
#' @param annotation A processed GRangesList or CompressedGRangesList object 
#' as returned by \code{process_input()}. 
#' @param inter_dir Output directory.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param ... Any additional arguments to
#' `mcscanx`.
#'
#' @return Paths to .collinearity files.
#' @noRd
interspecies_synteny <- function(blast_inter = NULL, annotation = NULL,
                                 inter_dir = NULL, anchors = 5, 
                                 max_gaps = 25, verbose = FALSE, ...) {
    species <- names(annotation)
    unique_comp <- combn(species, 2, simplify = FALSE)
    
    # Create BLAST tables and annotation data frames to export
    minter <- lapply(unique_comp, function(x) {
        n1 <- x[1]
        n2 <- x[2]
        amerged <- suppressWarnings(c(annotation[[n1]], annotation[[n2]]))
        amerged <- as.data.frame(amerged)
        amerged <- amerged[, c("seqnames", "gene", "start", "end")]
        blast1 <- paste0(n1, "_", n2)
        blast2 <- paste0(n2, "_", n1)
        bmerged <- rbind(blast_inter[[blast1]], blast_inter[[blast2]])
        
        res_list <- list(blast_table = bmerged, gff = amerged)
        return(res_list)
    })
    names(minter) <- unlist(lapply(unique_comp, function(x) {
        return(paste0(x[1], "_", x[2]))
    }))
    # Detect synteny
    intersyn <- unlist(lapply(seq_along(minter), function(x) {
        sp <- names(minter)[x]
        blast_file <- file.path(inter_dir, paste0(sp, ".blast"))
        write.table(minter[[x]]$blast_table, file = blast_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        gff_file <- file.path(inter_dir, paste0(sp, ".gff"))
        write.table(minter[[x]]$gff, file = gff_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        # Detect synteny
        input <- file.path(inter_dir, sp)
        if(verbose) { verbose <- "" }
        #syn_args <- c("-a -b 2 ", input, "-s ", anchors, "-m ", max_gaps)
        #syn <- system2("MCScanX", args = syn_args, stdout = verbose)
        rcpp_mcscanx_file(blast_file=blast_file, gff_file=gff_file, prefix=sp,
            outdir=inter_dir, match_size=5, max_gaps=max_gaps, is_pairwise=TRUE,
            in_synteny=2, ...)
        
        # Delete intermediate files
        unlink(c(blast_file, gff_file))
        synf <- file.path(inter_dir, paste0(sp, ".collinearity"))
        return(synf)
    }))
    intersyn <- intersyn[!is.null(intersyn)]
    return(intersyn)
}
