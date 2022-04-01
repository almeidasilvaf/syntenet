
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

