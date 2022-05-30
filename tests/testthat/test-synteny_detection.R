
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


test_that("intraspecies_synteny() detects intraspecies synteny", {
    blast_intra <- blast_list[1]
    annot <- as.data.frame(annotation[[1]])
    annot <- annot[, c("seqnames", "gene_id", "start", "end")]
    annot_list <- list(Olucimarinus = annot)
    intrasyn <- intraspecies_synteny(blast_intra, annot_list = annot_list)
    
    expect_equal(class(intrasyn), "character")
    expect_equal(length(intrasyn), 1)
})

test_that("interspecies_synteny() detects interspecies synteny", {
    data(blast_list)
    data(annotation)
    blast_inter <- blast_list[2]
    annot_list <- lapply(annotation, function(x) {
        x$gene <- x$gene_id
        return(x)
    })
    intersyn <- interspecies_synteny(blast_inter, annot_list = annot_list)
    
    expect_equal(class(intersyn), "character")
    expect_equal(length(intersyn), 1)
})

test_that("parse_collinearity() reads MCScanX output", {
    file <- system.file("extdata", "Olu.collinearity", package = "syntenet")
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
