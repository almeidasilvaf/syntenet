
#----Load data------------------------------------------------------------------
data(proteomes)
data(annotation)

#----Start tests----------------------------------------------------------------
test_that("check_input() returns TRUE if everything is OK", {
    c <- check_input(proteomes, annotation)
    expect_equal(c, TRUE)
    expect_equal(class(c), "logical")
})

test_that("process_input() returns clean seq and annotation objects", {
    clean <- process_input(proteomes, annotation)
    expect_equal(class(clean), "list")
    expect_false(identical(clean$seq, proteomes))
    expect_false(identical(clean$annotation, annotation))
})
