
#----Setup----------------------------------------------------------------------
data(proteomes)
data(annotation)

#----Start tests----------------------------------------------------------------
test_that("run_diamond() runs DIAMOND on the background", {
    seq <- process_input(proteomes, annotation)$seq[1:2]
    if(diamond_is_installed()) {
        diamond <- run_diamond(seq)
    } else {
        diamond <- list(
            A = data.frame(),
            B = data.frame(),
            C = data.frame(),
            D = data.frame()
        )
    }
    expect_equal(class(diamond), "list")
    expect_equal(length(diamond), length(seq) ^ 2)
    expect_equal(class(diamond[[1]]), "data.frame")
})


