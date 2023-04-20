
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
