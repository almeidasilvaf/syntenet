
#----Setup----------------------------------------------------------------------
data(proteomes)
data(annotation)
data(blast_list)
data(edges)

#----Start tests----------------------------------------------------------------
test_that("get_comp() returns a data frame of comparisons", {
    
    species <- c("sp1", "sp2")
    comp_df <- data.frame(sp1 = "sp1", sp2 = "sp2")
    t1 <- get_comp(species, compare = "all")
    t2 <- get_comp(species, compare = "intraspecies")
    t3 <- get_comp(species, compare = "interspecies")
    t4 <- get_comp(species, compare = comp_df)
    
    expect_error(
        get_comp(species, compare = cbind(comp_df, comp_df))
    )
    
    expect_error(
        get_comp(species, compare = data.frame(sp1 = "sp1", sp2 = "sp3"))
    )
    
    expect_error(
        get_comp(species, compare = "error")
    )
    
    expect_equal(nrow(t1), 4)
    expect_equal(nrow(t2), 2)
    expect_equal(nrow(t3), 2)
    expect_equal(nrow(t4), 1)
    expect_true(all(names(t1) %in% c("sp1", "sp2")))
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


test_that("intraspecies_synteny() detects intraspecies synteny", {
    
    pdata <- process_input(proteomes, annotation)
    pannotation <- pdata$annotation[1]
    blast_intra <- blast_list[1]

    # Detect intraspecies synteny
    intrasyn <- intraspecies_synteny(blast_intra, pannotation)

    expect_equal(class(intrasyn), "character")
    expect_equal(length(intrasyn), 1)
})

test_that("interspecies_synteny() detects interspecies synteny", {
    
    # Get DIAMOND and processed annotation lists
    blast_inter <- blast_list[2]
    pannotation <- process_input(proteomes, annotation)$annotation

    # Detect interspecies synteny
    intersyn <- interspecies_synteny(blast_inter, pannotation)
    
    expect_equal(class(intersyn), "character")
    expect_equal(length(intersyn), 1)
})

test_that("parse_collinearity() reads MCScanX output", {
    file <- system.file("extdata", "Scerevisiae.collinearity", package = "syntenet")
    net <- parse_collinearity(file)
    
    expect_equal(class(net), "data.frame")
    expect_equal(ncol(net), 2)
})


test_that("infer_syntenet() on example set find anchors", {
    processed <- process_input(proteomes, annotation) 
    annotation <- processed$annotation
    net <- infer_syntenet(blast_list, annotation)
    
    expect_equal(net, edges)
})

test_that("infer_syntenet() on example set find anchors verbose", {
    processed <- process_input(proteomes, annotation) 
    seq <- processed$seq
    annotation <- processed$annotation
    net <- infer_syntenet(blast_list, annotation, verbose = TRUE)
    
    expect_equal(net, edges)
})

test_that("infer_syntenet() on example set write html", {
    processed <- process_input(proteomes, annotation) 
    seq <- processed$seq
    annotation <- processed$annotation
    net <- infer_syntenet(blast_list, annotation, is_pairwise = FALSE)
    
    expect_equal(net, edges)
})

test_that("infer_syntenet() on example set write html verbose", {
    processed <- process_input(proteomes, annotation) 
    seq <- processed$seq
    annotation <- processed$annotation
    net <- infer_syntenet(blast_list, annotation, is_pairwise = FALSE,
        verbose = TRUE)
    
    expect_equal(net, edges)
})
