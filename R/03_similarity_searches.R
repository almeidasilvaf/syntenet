

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
#' @param compare Character scalar indicating which comparisons
#' should be made when running DIAMOND. 
#' Possible modes are "all" (all-vs-all comparisons), 
#' "intraspecies" (intraspecies comparisons only), or
#' "interspecies" (interspecies comparisons only). Alternatively, users can
#' pass a 2-column data frame as input with the names of species to be 
#' compared.
#' @param ... Any additional arguments to `diamond blastp`.
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
                        outdir = tempdir(), threads = NULL, 
                        compare = "all", ...) {
    
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
    comb_df <- get_comp(names(seq), compare = compare)
    diamond_blastp <- lapply(seq_len(nrow(comb_df)), function(x) {
        query <- paste0(file.path(seqdir, comb_df[x, 1]), ".fasta")
        db <- file.path(dbdir, comb_df[x, 2])
        outfile <- paste(comb_df[x, 1], comb_df[x, 2], sep = "_")
        outfile <- paste0(file.path(resdir, outfile), sep =  ".tsv")
        
        bargs <- c("blastp -q", query, "-d", db, "-o", outfile, threads,
                   "--max-hsps 1 -k", top_hits, "--quiet", ...)
        run_diamond <- system2("diamond", args = bargs)
    })
    
    final_list <- read_diamond(resdir)
    return(final_list)
}


#' Wrapper to run last from an R session
#'
#' @param seq A processed list of AAStringSet objects 
#' as returned by \code{process_input()}.
#' @param verbose Logical indicating if progress messages should be printed.
#' Default: FALSE. 
#' @param outdir Output directory for last results. By default, output files
#' are saved to a temporary directory.
#' @param threads Number of threads to use. Default: 1.
#' @param compare Character scalar indicating which comparisons
#' should be made when running last. 
#' Possible modes are "all" (all-vs-all comparisons), 
#' "intraspecies" (intraspecies comparisons only), or
#' "interspecies" (interspecies comparisons only). Alternatively, users can
#' pass a 2-column data frame as input with the names of species to be 
#' compared.
#' @param lastD last option D: query letters per random alignment. Default: 1e6.
#' @param ... Any additional arguments to `lastal`.
#'
#' @return A list of data frames containing last's tabular output
#' for each pairwise combination of species. For n species, the list length
#' will be \eqn{n^2}.
#' @export
#' @rdname run_last
#' @importFrom Biostrings writeXStringSet
#' @importFrom utils read.csv
#' @examples 
#' data(proteomes)
#' data(annotation)
#' seq <- process_input(proteomes, annotation)$seq[1:2]
#' if(last_is_installed()) {
#'     last_results <- run_last(seq)
#' }
run_last <- function(seq = NULL, verbose = FALSE, 
                        outdir = tempdir(), threads = 1, 
                        compare = "all", lastD=1e6, ...) {
    
    valid <- valid_seq(seq)
    check_last <- last_is_installed()
    if(!check_last) {
        stop("Unable to find last in PATH.")
    }
    
    # 1. Make dbs for each species
    if(verbose) { message("1. Creating database for each species...\n") }
    dbdir <- file.path(outdir, "last", "dbs")
    seqdir <- file.path(outdir, "last", "seqs")
    resdir <- file.path(outdir, "last", "results")
    if(!dir.exists(dbdir)) { dir.create(dbdir, recursive = TRUE) }
    if(!dir.exists(seqdir)) { dir.create(seqdir, recursive = TRUE) }
    if(!dir.exists(resdir)) { dir.create(resdir, recursive = TRUE) }
    
    makedb <- lapply(seq_along(seq), function(x) {
        seqfile <- paste0(file.path(seqdir, names(seq)[x]), ".fasta")
        dbfile <- file.path(dbdir, names(seq)[x])
        Biostrings::writeXStringSet(seq[[x]], filepath = seqfile)
        dbargs <- c("-p", "-cR01", "-P", threads, dbfile, seqfile)
        rundb <- system2("lastdb", args = dbargs)
    })
    
    # 2. Pairwise BLASTp-like search with last
    if(verbose) { message("2. Running pairwise last searches...\n")}
    comb_df <- get_comp(names(seq), compare = compare)
    last_blastp <- lapply(seq_len(nrow(comb_df)), function(x) {
        query <- paste0(file.path(seqdir, comb_df[x, 1]), ".fasta")
        db <- file.path(dbdir, comb_df[x, 2])
        outfile <- paste(comb_df[x, 1], comb_df[x, 2], sep = "_")
        outfile <- paste0(file.path(resdir, outfile), sep =  ".tsv")
        
        bargs <- c("-f", "BlastTab+", "-P", threads, "-D", lastD, db, query, ..., ">", outfile)
        run_last <- system2("lastal", args = bargs)
    })
    
    final_list <- read_last(resdir)
    return(final_list)
}


#' Export processed sequences as FASTA files
#' 
#' @param seq A processed list of AAStringSet objects 
#' as returned by \code{process_input()}.
#' @param outdir Path to output directory where FASTA files will be stored.
#' 
#' @return Path to exported FASTA files.
#' 
#' @rdname export_sequences
#' @importFrom Biostrings writeXStringSet
#' @export
#' @examples
#' # Load data
#' data(proteomes)
#' data(annotation)
#' 
#' # Process data
#' pdata <- process_input(proteomes, annotation)
#' 
#' # Export data
#' outdir <- file.path(tempdir(), "example_test")
#' export_sequences(pdata$seq, outdir)
export_sequences <- function(seq = NULL, outdir = tempdir()) {
    
    if(!dir.exists(outdir)) { dir.create(outdir, recursive = TRUE) }
    
    # Iterate through sequences and export them
    paths <- unlist(lapply(seq_along(seq), function(x) {
        
        filename <- paste0(file.path(outdir, names(seq)[x]), ".fasta")
        w <- Biostrings::writeXStringSet(
            seq[[x]], 
            filepath = filename,
            compress = FALSE
        )
        return(filename)
    }))
    
    return(paths)
}


#' Read DIAMOND/BLAST tables as a list of data frames
#'
#' @param diamond_dir Path to directory containing the tabular output
#' of DIAMOND or similar programs (e.g., BLAST).
#'
#' @return A list of data frames with the tabular DIAMOND output.
#' 
#' @rdname read_diamond
#' @export
#' @importFrom utils read.csv
#' @examples
#' # Path to output directory
#' diamond_dir <- system.file("extdata", package = "syntenet")
#' 
#' # Read output
#' l <- read_diamond(diamond_dir)
read_diamond <- function(diamond_dir = NULL) {
    
    if(!dir.exists(diamond_dir) | is.null(diamond_dir)) {
        stop("Could not find the directory specified in 'diamond_dir'.")
    }
    
    # Read files as a list of data frames
    result_files <- list.files(diamond_dir, pattern = ".tsv", full.names = TRUE)
    final_list <- lapply(result_files, function(x) {
        df <- read.csv(x, header = FALSE, sep = "\t")
        names(df) <- c(
            "query", "db", "perc_identity", "length", "mismatches", 
            "gap_open", "qstart", "qend", "tstart", "tend",
            "evalue", "bitscore"
        )
        return(df)
    })
    names(final_list) <- gsub("\\.tsv", "", basename(result_files))
    
    return(final_list)
}


#' Read last/BLAST tables as a list of data frames
#'
#' @param last_dir Path to directory containing the tabular output
#' of last or similar programs (e.g., BLAST).
#'
#' @return A list of data frames with the tabular last output.
#' 
#' @rdname read_last
#' @export
#' @importFrom utils read.csv
#' @examples
#' # Path to output directory
#' last_dir <- system.file("extdata", package = "syntenet")
#' 
#' # Read output
#' l <- read_last(last_dir)
read_last <- function(last_dir = NULL) {
    
    if(!dir.exists(last_dir) | is.null(last_dir)) {
        stop("Could not find the directory specified in 'last_dir'.")
    }
    
    # Read files as a list of data frames
    result_files <- list.files(last_dir, pattern = ".tsv", full.names = TRUE)
    final_list <- lapply(result_files, function(x) {
        df <- read.csv(x, header = FALSE, sep = "\t")
        names(df) <- c(
            "query", "db", "perc_identity", "length", "mismatches", 
            "gap_open", "qstart", "qend", "tstart", "tend",
            "evalue", "bitscore", "qlen", "tlen", "rawscore"
        )
        return(df)
    })
    names(final_list) <- gsub("\\.tsv", "", basename(result_files))
    
    return(final_list)
}

