
#' Detect intraspecies synteny
#'
#' @param blast_intra A list of BLAST data frames for intraspecies comparisons.
#' @param intra_dir Output directory.
#' @param annot_list A list of annotation data frames.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param is_pairwise Logical indicating if only pairwise blocks should be 
#' reported. Default: TRUE.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param ... Any additional arguments to
#' `mcscanx`.
#'
#' @return Paths to .collinearity files.
#' @rdname intraspecies_synteny
#' @export
#' @examples 
#' data(blast_list)
#' data(annotation)
#' blast_intra <- blast_list[1]
#' annot <- as.data.frame(annotation[[1]])
#' annot <- annot[, c("seqnames", "gene_id", "start", "end")]
#' annot_list <- list(Olucimarinus = annot)
#' intrasyn <- intraspecies_synteny(blast_intra, annot_list = annot_list)
intraspecies_synteny <- function(blast_intra = NULL, 
                                 intra_dir = file.path(tempdir(), "intra"),
                                 annot_list = NULL, anchors = 5, 
                                 max_gaps = 25, is_pairwise = TRUE,
                                 verbose = FALSE, ...) {

    if(!dir.exists(intra_dir)) { dir.create(intra_dir, recursive = TRUE) }
    
    intrasyn <- unlist(lapply(blast_intra, function(x) {
        sp <- gsub("_.*", "", x[1, 1])
        
        # Write files
        blast_file <- file.path(intra_dir, paste0(sp, ".blast"))
        write.table(x, file = blast_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        gff_file <- file.path(intra_dir, paste0(sp, ".gff"))
        write.table(annot_list[[sp]], file = gff_file, quote = FALSE,
                    row.names = FALSE, col.names = FALSE, sep = "\t")
        
        # Detect synteny
        input <- file.path(intra_dir, sp)
        rcpp_mcscanx_file(
            blast_file = blast_file, gff_file = gff_file, prefix = sp,
            outdir = intra_dir, match_size = anchors, max_gaps = max_gaps,
            is_pairwise = is_pairwise, in_synteny = 1, verbose = verbose, ...
        )
        
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
#' @param annot_list A processed GRangesList or CompressedGRangesList object 
#' as returned by \code{process_input()}. 
#' @param inter_dir Output directory.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param is_pairwise specify if only pairwise blocks should be reported
#' Default: TRUE.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param ... Any additional arguments to
#' `mcscanx`.
#'
#' @return Paths to .collinearity files.
#' 
#' @importFrom utils combn
#' @rdname interspecies_synteny
#' @export
#' @examples 
#' data(blast_list)
#' data(annotation)
#' blast_inter <- blast_list[2]
#' annot_list <- lapply(annotation, function(x) {
#'     x$gene <- x$gene_id
#'     return(x)
#' })
#' intersyn <- interspecies_synteny(blast_inter, annot_list = annot_list)
interspecies_synteny <- function(blast_inter = NULL, annot_list = NULL,
                                 inter_dir = file.path(tempdir(), "inter"), 
                                 anchors = 5, 
                                 max_gaps = 25, is_pairwise = TRUE,
                                 verbose = FALSE, ...) {
    
    if(!dir.exists(inter_dir)) { dir.create(inter_dir, recursive = TRUE) }
    species <- names(annot_list)
    unique_comp <- combn(species, 2, simplify = FALSE)
    
    # Create BLAST tables and annotation data frames to export
    minter <- lapply(unique_comp, function(x) {
        n1 <- x[1]
        n2 <- x[2]
        amerged <- suppressWarnings(c(annot_list[[n1]], annot_list[[n2]]))
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
        rcpp_mcscanx_file(
            blast_file = blast_file, gff_file = gff_file, prefix = sp,
            outdir = inter_dir, match_size = anchors, max_gaps = max_gaps,
            is_pairwise = is_pairwise, in_synteny = 2, verbose = verbose, ...
        )
        
        # Delete intermediate files
        unlink(c(blast_file, gff_file))
        synf <- file.path(inter_dir, paste0(sp, ".collinearity"))
        return(synf)
    }))
    intersyn <- intersyn[!is.null(intersyn)]
    return(intersyn)
}
