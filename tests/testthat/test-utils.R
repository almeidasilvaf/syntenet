
#----Load data------------------------------------------------------------------
data(proteomes)
data(annotation)

#----Start tests----------------------------------------------------------------
test_that("check_list_names() checks if names match", {
    check <- check_list_names(proteomes, annotation)
    expect_equal(check, TRUE)
})


test_that("check_ngenes() checks if number of seqs is <= gene count", {
    check <- check_ngenes(proteomes, annotation)
    expect_equal(check, TRUE)
})

test_that("check_gene_names() checks if sequence names match gene names", {
    check <- check_gene_names(proteomes, annotation)
    expect_equal(check, TRUE)
})

test_that("species_id_length() returns the best string length for IDs", {
    len <- species_id_length(proteomes)
    expect_equal(class(len), "numeric")
})

test_that("is_valid() returns a logical scalar", {
    valid <- is_valid("echo", "Hello world")
    expect_equal(class(valid), "logical")
})

test_that("diamond_is_installed() returns a logical scalar", {
    diamond <- diamond_is_installed()
    expect_equal(class(diamond), "logical")
})
