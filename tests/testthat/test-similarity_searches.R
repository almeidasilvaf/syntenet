
# Load data
data(proteomes)
data(annotation)
data(blast_list)

# Start tests
test_that("make_bidirectional() creates bidirectional comparisons", {
    spp_outgroup <- data.frame(species = "spA", outgroup = "spB")
    comp_bi <- make_bidirectional(spp_outgroup)
    
    n1 <- as.integer(nrow(spp_outgroup))
    n2 <- as.integer(nrow(comp_bi))
    expect_true(n2 == 2L*n1)
})


test_that("collapse_bidirectional_hits() collapses hits", {
    
    blast_inter <- blast_list[c(2,3)]
    compare <- data.frame(species = "Olucimarinus", outgroup = "OspRCC809")
    chits <- collapse_bidirectional_hits(blast_inter, compare)
    
    expect_equal(length(chits), 1)
})


test_that("run_diamond() runs DIAMOND on the background", {
    
    seq <- process_input(proteomes, annotation)$seq[1:2]
    if(diamond_is_installed()) {
        diamond <- run_diamond(seq, verbose = TRUE, threads = 1)
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


test_that("export_sequences() exports processed sequences to FASTA files", {
    
    # Process data
    pdata <- process_input(proteomes, annotation)

    # Export data
    outdir <- file.path(tempdir(), "example_testthat")
    p <- export_sequences(pdata$seq, outdir)
    
    expect_equal(class(p), "character")
    
})


test_that("read_diamond() creates a list of data frames", {
    
    diamond_dir <- system.file("extdata", package = "syntenet")

    # Read output
    l <- read_diamond(diamond_dir)
    
    expect_error(read_diamond())
    expect_equal(class(l), "list")
    expect_equal(class(l[[1]]), "data.frame")
})

