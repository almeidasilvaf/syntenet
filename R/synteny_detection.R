
#' Wrapper to run DIAMOND from an R session
#'
#' @param seq A processed list of AAStringSet objects 
#' as returned by \code{process_input()}.
#' @param top_hits Number of top hits to keep in DIAMOND search. Default: 5.
#' @param verbose Logical indicating if progress messages should be printed.
#' Default: FALSE. 
#' @param outdir Output directory for DIAMOND results. By default, output files
#' are saved to a temporary directory.
#' @param threads Number of threads to use. Default: let DIAMOND auto-detect
#' and use all available virtual cores on the machine.
#' @param ... Any additional arguments to
#' `diamond blastp`.
#'
#' @return A list of data frames containing DIAMOND's tabular output
#' for each pairwise combination of species. For n species, the list length
#' will be \eqn{n^2}.
#' @export
#' @rdname run_diamond
#' @importFrom Biostrings writeXStringSet
#' @importFrom utils read.csv
#' @examples 
#' data(proteomes)
#' data(annotation)
#' seq <- process_input(proteomes, annotation)$seq[1:2]
#' if(diamond_is_installed()) {
#'     diamond_results <- run_diamond(seq)
#' }
run_diamond <- function(seq = NULL, top_hits = 5, verbose = FALSE, 
                        outdir = tempdir(), threads = NULL, ...) {
    
    valid <- valid_seq(seq)
    check_diamond <- diamond_is_installed()
    if(!check_diamond) {
        stop("Unable to find DIAMOND in PATH.")
    }
    
    # 1. Make dbs for each species
    if(verbose) { message("1. Creating database for each species...\n") }
    dbdir <- file.path(outdir, "diamond", "dbs")
    seqdir <- file.path(outdir, "diamond", "seqs")
    resdir <- file.path(outdir, "diamond", "results")
    if(!dir.exists(dbdir)) { dir.create(dbdir, recursive = TRUE) }
    if(!dir.exists(seqdir)) { dir.create(seqdir, recursive = TRUE) }
    if(!dir.exists(resdir)) { dir.create(resdir, recursive = TRUE) }
    
    makedb <- lapply(seq_along(seq), function(x) {
        seqfile <- paste0(file.path(seqdir, names(seq)[x]), ".fasta")
        dbfile <- file.path(dbdir, names(seq)[x])
        Biostrings::writeXStringSet(seq[[x]], filepath = seqfile)
        dbargs <- c("makedb --in ", seqfile, "-d ", dbfile, "--quiet")
        rundb <- system2("diamond", args = dbargs)
    })
    
    if(!is.null(threads)) { threads <- paste0("-p ", threads) }
    # 2. Pairwise BLASTp-like search with DIAMOND
    if(verbose) { message("2. Running pairwise DIAMOND searches...\n")}
    comb_df <- expand.grid(names(seq), names(seq), stringsAsFactors = FALSE)
    diamond_blastp <- lapply(seq_len(nrow(comb_df)), function(x) {
        query <- paste0(file.path(seqdir, comb_df[x, 1]), ".fasta")
        db <- file.path(dbdir, comb_df[x, 2])
        outfile <- paste(comb_df[x, 1], comb_df[x, 2], sep = "_")
        outfile <- paste0(file.path(resdir, outfile), sep =  ".tsv")
        
        bargs <- c("blastp -q", query, "-d", db, "-o", outfile, threads,
                   "--max-hsps 1 -k", top_hits, "--quiet", ...)
        run_diamond <- system2("diamond", args = bargs)
    })
    
    result_files <- list.files(resdir, pattern = ".tsv", full.names = TRUE)
    final_list <- lapply(result_files, function(x) {
        df <- utils::read.csv(x, header = FALSE, sep = "\t")
        names(df) <- c("query", "db", "perc_identity", "length", "mismatches", 
                       "gap_open", "qstart", "qend", "tstart", "tend",
                       "evalue", "bitscore")
        return(df)
    })
    names(final_list) <- gsub("\\.tsv", "", basename(result_files))
    return(final_list)
}



#' Parse .collinearity files into an edge list
#'
#' @param collinearity_paths Character vector of paths to .collinearity files.
#'
#' @return A 2-column data frame with each gene from the anchor pair in each
#' column.
#'
#' @importFrom utils read.table
#' @export
#' @rdname parse_collinearity
#' @examples 
#' collinearity_paths <- system.file("extdata", "Olu.collinearity", 
#'                                   package = "syntenet")
#' net <- parse_collinearity(collinearity_paths)
parse_collinearity <- function(collinearity_paths = NULL) {
    
    fname <- gsub("\\.collinearity", "", basename(collinearity_paths))
    names(collinearity_paths) <- fname
    
    netdb <- lapply(seq_along(collinearity_paths), function(x) {
        lines <- readLines(collinearity_paths[x])
        nlines <- length(lines[!startsWith(lines, "#")])
        
        df <- NULL
        if(nlines > 0) {
            df <- read.table(collinearity_paths[x], sep = "\t", comment.char = "#")
            df <- df[, c(2,3)]
            names(df) <- c("Anchor1", "Anchor2")
        }
        return(df)
    })
    netdb <- Reduce(rbind, netdb)
    return(netdb)
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
#' data(proteomes)
#' data(annotation)
#' data(blast_list)
#' processed <- process_input(proteomes, annotation) 
#' seq <- processed$seq
#' annotation <- processed$annotation
#' net <- infer_syntenet(blast_list, annotation)
infer_syntenet <- function(blast_list = NULL, annotation = NULL,
                           outdir = tempdir(),
                           anchors = 5, max_gaps = 25,
                           is_pairwise = TRUE, verbose = FALSE, ...) {

    valid <- valid_blast(blast_list) & valid_annot(annotation)
    annot_dfs <- lapply(annotation, function(x) {
        return(as.data.frame(x)[, c("seqnames", "gene", "start", "end")])
    })
    names(annot_dfs) <- unlist(lapply(annotation, function(x) {
        return(gsub("_.*", "", x$gene[1]))
    }))

    # Separate intra from interspecies
    species <- names(annotation)
    equal_comp <- vapply(species, function(x) paste0(x, "_", x), character(1))
    idx_equal <- which(names(blast_list) %in% equal_comp)
    
    #---- 1) Intraspecies synteny detection------------------------------------
    intra_dir <- file.path(outdir, "intraspecies_synteny")
    if(!dir.exists(intra_dir)) { dir.create(intra_dir, recursive = TRUE) }
    blast_intra <- blast_list[idx_equal]
    
    intraspecies <- intraspecies_synteny(blast_intra, intra_dir, annot_dfs,
                                         anchors, max_gaps, is_pairwise,
                                         verbose, ...)
    
    #---- 2) Interspecies synteny detection------------------------------------
    inter_dir <- file.path(outdir, "interspecies_synteny")
    if(!dir.exists(inter_dir)) { dir.create(inter_dir, recursive = TRUE) }
    blast_inter <- blast_list[-idx_equal]
    
    interspecies <- interspecies_synteny(blast_inter, annotation, inter_dir,
                                         anchors, max_gaps, is_pairwise,
                                         verbose, ...)
    
    # Create edge list
    syn_files <- c(intraspecies, interspecies)
    edges <- parse_collinearity(syn_files)
    return(edges)
}
