
#----Setup----------------------------------------------------------------------
data(proteomes)
data(annotation)
data(blast_list)
data(edges)

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

test_that("parse_collinearity() reads MCScanX output", {
    file <- system.file("extdata", "Olu.collinearity", package = "syntenet")
    net <- parse_collinearity(file)
    
    expect_equal(class(net), "data.frame")
    expect_equal(ncol(net), 2)
})

test_that("infer_syntenet() on example set find anchors", {
    processed <- process_input(proteomes, annotation) 
    seq <- processed$seq
    annotation <- processed$annotation
    net <- infer_syntenet(blast_list, annotation)
    
    expect_equal(net, edges)
})
