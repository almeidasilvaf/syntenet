
#' Detect intraspecies synteny
#'
#' @param blast_intra A list of BLAST/DIAMOND data frames for 
#' intraspecies comparisons as returned by \code{run_diamond()}.
#' @param annotation A processed GRangesList or CompressedGRangesList object 
#' as returned by \code{process_input()}. 
#' @param intra_dir Path to output directory where .collinearity files will
#' be stored.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param is_pairwise Logical indicating if only pairwise blocks should be 
#' reported. Default: TRUE.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param bp_param BiocParallel back-end to be used. 
#' Default: \code{BiocParallel::SerialParam()}.
#' @param ... Any additional arguments to the `MCScanX` algorithm. For a
#' complete list of all available options, see the man page 
#' of \code{rcpp_mcscanx_file()}.
#'
#' @return Paths to .collinearity files.
#' @importFrom methods is
#' @importFrom BiocParallel SerialParam bplapply
#' @rdname intraspecies_synteny
#' @export
#' @examples 
#' # Load data
#' data(scerevisiae_annot)
#' data(scerevisiae_diamond)
#' 
#' # Detect intragenome synteny
#' intra_syn <- intraspecies_synteny(
#'     scerevisiae_diamond, scerevisiae_annot
#' )
#' 
intraspecies_synteny <- function(
        blast_intra = NULL, 
        annotation = NULL, 
        intra_dir = file.path(tempdir(), "intra"),
        anchors = 5, max_gaps = 25, is_pairwise = TRUE,
        verbose = FALSE, 
        bp_param = BiocParallel::SerialParam(),
        ...
) {
    
    valid <- valid_blast(blast_intra) & valid_annot(annotation)
    if(!dir.exists(intra_dir)) { dir.create(intra_dir, recursive = TRUE) }
    
    # Get species ID length
    example_gene <- as.character(blast_intra[[1]][1, 1])
    species_id_length <- nchar(gsub("_.*", "", example_gene))
    
    intrasyn <- BiocParallel::bplapply(seq_along(blast_intra), function(x) {
        
        # Get species name
        blast_comp <- names(blast_intra)[x]
        all_species <- names(annotation)
        counts <- unlist(lapply(all_species, function(y) { 
            match1 <- grepl(paste0(y, "_"), blast_comp) 
            match2 <- grepl(paste0(y, "$"), blast_comp)
            return(match1 + match2)
        }))
        species <- all_species[counts == 2]
        
        # Write files
        ## BLAST output
        blast_file <- file.path(intra_dir, paste0(species, ".blast"))
        write.table(
            blast_intra[[x]], file = blast_file, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t"
        )
        
        ## Annotation (as a 4-column table)
        annot_df <- as.data.frame(annotation[[species]])
        annot_df <- annot_df[, c("seqnames", "gene", "start", "end")]
        
        gff_file <- file.path(intra_dir, paste0(species, ".gff"))
        write.table(
            annot_df, file = gff_file, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t"
        )
        
        # Detect synteny
        input <- file.path(intra_dir, species)
        rcpp_mcscanx_file(
            blast_file = blast_file, gff_file = gff_file, prefix = species,
            outdir = intra_dir, match_size = anchors, max_gaps = max_gaps,
            is_pairwise = is_pairwise, in_synteny = 1, verbose = verbose, 
            species_id_length = species_id_length, ...
        )
        
        # Delete intermediate files
        unlink(c(blast_file, gff_file))
        
        synf <- file.path(intra_dir, paste0(species, ".collinearity"))
        
        return(synf)
    }, BPPARAM = bp_param)
    intrasyn <- unlist(intrasyn)
    intrasyn <- intrasyn[!is.null(intrasyn)]
    
    return(intrasyn)
}


#' Detect interspecies synteny
#'
#' @param blast_inter A list of BLAST/DIAMOND data frames for 
#' interspecies comparisons as returned by \code{run_diamond()}.
#' @param annotation A processed GRangesList or CompressedGRangesList object 
#' as returned by \code{process_input()}.
#' @param inter_dir Path to output directory where .collinearity files will
#' be stored.
#' @param anchors Numeric indicating the minimum required number of genes
#' to call a syntenic block. Default: 5.
#' @param max_gaps Numeric indicating the number of upstream and downstream
#' genes to search for anchors. Default: 25.
#' @param is_pairwise specify if only pairwise blocks should be reported
#' Default: TRUE.
#' @param verbose Logical indicating if log messages should be printed on
#' screen. Default: FALSE.
#' @param bp_param BiocParallel back-end to be used. 
#' Default: \code{BiocParallel::SerialParam()}.
#' @param ... Any additional arguments to the `MCScanX` algorithm. For a
#' complete list of all available options, see the man page 
#' of \code{rcpp_mcscanx_file()}.
#'
#' @return Paths to .collinearity files.
#' 
#' @importFrom utils combn
#' @importFrom methods is
#' @importFrom BiocParallel bplapply SerialParam
#' @rdname interspecies_synteny
#' @export
#' @examples 
#' # Load data
#' data(proteomes)
#' data(blast_list)
#' data(annotation)
#' 
#' # Get DIAMOND and processed annotation lists
#' blast_inter <- blast_list[2]
#' annotation <- process_input(proteomes, annotation)$annotation
#'
#' # Detect interspecies synteny
#' intersyn <- interspecies_synteny(blast_inter, annotation)
interspecies_synteny <- function(
        blast_inter = NULL, 
        annotation = NULL,
        inter_dir = file.path(tempdir(), "inter"), 
        anchors = 5, max_gaps = 25, is_pairwise = TRUE,
        verbose = FALSE, 
        bp_param = BiocParallel::SerialParam(),
        ...
) {
    
    valid <- valid_blast(blast_inter) & valid_annot(annotation)
    if(!dir.exists(inter_dir)) { dir.create(inter_dir, recursive = TRUE) }
    species <- names(annotation)
    unique_comp <- combn(species, 2, simplify = FALSE)
    
    # Create a list with BLAST tables and annotation data frames to export
    minter <- lapply(unique_comp, function(x) {
        n1 <- x[1]
        n2 <- x[2]
        
        ## rbind GRanges into a single annotation data frame
        amerged <- suppressWarnings(c(annotation[[n1]], annotation[[n2]]))
        amerged <- as.data.frame(amerged)
        amerged <- amerged[, c("seqnames", "gene", "start", "end")]
        
        ## rbind BLAST tables into a single data frame
        blast1 <- paste0(n1, "_", n2)
        blast2 <- paste0(n2, "_", n1)
        bmerged <- rbind(blast_inter[[blast1]], blast_inter[[blast2]])
        
        res_list <- list(blast_table = bmerged, gff = amerged)
        return(res_list)
    })
    names(minter) <- unlist(lapply(unique_comp, function(x) {
        return(paste0(x[1], "_", x[2]))
    }))
    
    # Get species ID length
    example_chr <- as.character(minter[[1]]$gff[1, 1])
    species_id_length <- nchar(gsub("_.*", "", example_chr))
    
    # Detect synteny
    intersyn <- BiocParallel::bplapply(seq_along(minter), function(x) {
        
        sp <- names(minter)[x]
        blast_file <- file.path(inter_dir, paste0(sp, ".blast"))
        write.table(
            minter[[x]]$blast_table, file = blast_file, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t"
        )
        
        gff_file <- file.path(inter_dir, paste0(sp, ".gff"))
        write.table(
            minter[[x]]$gff, file = gff_file, quote = FALSE,
            row.names = FALSE, col.names = FALSE, sep = "\t"
        )
        
        # Detect synteny
        input <- file.path(inter_dir, sp)
        rcpp_mcscanx_file(
            blast_file = blast_file, gff_file = gff_file, prefix = sp,
            outdir = inter_dir, match_size = anchors, max_gaps = max_gaps,
            is_pairwise = is_pairwise, in_synteny = 2, verbose = verbose, 
            species_id_length = species_id_length, ...
        )
        
        # Delete intermediate files
        unlink(c(blast_file, gff_file))
        synf <- file.path(inter_dir, paste0(sp, ".collinearity"))
        return(synf)
    }, BPPARAM = bp_param)
    intersyn <- unlist(intersyn)
    intersyn <- intersyn[!is.null(intersyn)]
    
    return(intersyn)
}
