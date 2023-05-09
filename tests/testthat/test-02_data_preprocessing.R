
#----Load data------------------------------------------------------------------
data(proteomes)
data(annotation)

## Simulate a GRangesList object without "gene_id"
annotation2 <- lapply(annotation, function(x) {
    x$gene_id <- NULL
    return(x)
})

#----Start tests----------------------------------------------------------------
test_that("check_input() returns TRUE if everything is OK", {
    
    c <- check_input(proteomes, annotation)
    
    expect_equal(c, TRUE)
    expect_equal(class(c), "logical")
})

test_that("process_input() returns clean seq and annotation objects", {
    
    clean <- process_input(proteomes, annotation)
    expect_error(process_input(proteomes, annotation2))
    
    expect_equal(class(clean), "list")
    expect_false(identical(clean$seq, proteomes))
    expect_false(identical(clean$annotation, annotation))

})


test_that("collapse_protein_ids() works", {

    seq_path <- system.file(
        "extdata", "RefSeq_parsing_example", package = "syntenet"
    )
    seq <- fasta2AAStringSetlist(seq_path)
    annot <- gff2GRangesList(seq_path)
    names(seq$Aalosa) <- gsub(" .*", "", names(seq$Aalosa))
    cor_df <- as.data.frame(annot$Aalosa[annot$Aalosa$type == "CDS", ])
    cor_df <- cor_df[, c("Name", "gene")]
    protein2gene <- list(Aalosa = cor_df)
    new_seqs <- collapse_protein_ids(seq, protein2gene)
    
    # Start tests
    expect_true(is(new_seqs, "list"))
    expect_true(is(new_seqs[[1]], "AAStringSet"))
    
    expect_error(collapse_protein_ids(seq, protein2gene$Aalosa))
    
    seq2 <- seq
    names(seq2) <- "random_name"
    expect_error(collapse_protein_ids(seq2, protein2gene))
    
    pgene <- protein2gene
    pgene$Aalosa$Name <- pgene$Aalosa$gene
    expect_error(collapse_protein_ids(seq, pgene))
})