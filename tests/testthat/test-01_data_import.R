
#----Load data------------------------------------------------------------------
fasta_dir <- system.file("extdata", "sequences", package = "syntenet")
gff_dir <- system.file("extdata", "annotation", package = "syntenet")

#----Start tests----------------------------------------------------------------
test_that("gff2GRangesList() reads GFF/GTF files as a GRangesList object", {
    
    grangeslist <- gff2GRangesList(gff_dir)
    
    expect_true(grepl("GRangesList", class(grangeslist)))
})

test_that("fasta2AAStringSetlist() reads FASTA files as list of AAStringSet", {
    
    aastringsetlist <- fasta2AAStringSetlist(fasta_dir)
    
    expect_equal(class(aastringsetlist), "list")
    expect_true("AAStringSet" %in% class(aastringsetlist[[1]]))
})
