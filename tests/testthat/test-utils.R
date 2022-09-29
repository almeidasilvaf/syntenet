
#----Load data------------------------------------------------------------------
data(proteomes)
data(annotation)
data(blast_list)

## Simulate `proteomes` object without names
proteomes2 <- proteomes
names(proteomes2) <- NULL

## Simulate `proteomes` with names that do not match names in `annotation`
proteomes3 <- proteomes
names(proteomes3)[1] <- "different_name"

## Simulate `annotation` object with less genes than in `proteomes`
annotation2 <- annotation
annotation2$Olucimarinus <- annotation2$Olucimarinus[1:500]

## Simulate `annotation` without "Name" and "gene_id" columns
annotation3 <- annotation
annotation3$Olucimarinus$Name <- NULL
annotation3$Olucimarinus$gene_id <- NULL

## Simulate `annotation` with gene names that are different from `proteomes`
annotation4 <- annotation2
annotation4$Olucimarinus$gene_id[1:5] <- paste0("fakegene", 1:5) 

#----Start tests----------------------------------------------------------------
test_that("check_list_names() checks if names match", {
    
    check <- check_list_names(proteomes, annotation)

    expect_error(
        check_list_names(proteomes2, annotation)
    )
    
    expect_error(
        check_list_names(proteomes3, annotation)
    )
    
    expect_equal(check, TRUE)
})


test_that("check_ngenes() checks if number of seqs is <= gene count", {
    
    check <- check_ngenes(proteomes, annotation)
    
    expect_error(
        check_ngenes(proteomes, annotation2)
    )
    
    expect_equal(check, TRUE)
})


test_that("check_gene_names() checks if sequence names match gene names", {
    
    check <- check_gene_names(proteomes, annotation)
    
    expect_error(
        check_gene_names(proteomes, annotation3)
    )
    
    expect_error(
        check_gene_names(proteomes, annotation4)
    )
    
    expect_equal(check, TRUE)
})


test_that("species_id_length() returns the best string length for IDs", {
    len <- species_id_length(proteomes)
    expect_true(class(len) %in% c("numeric", "integer"))
})


test_that("is_valid() returns a logical scalar", {
    valid <- is_valid("echo", "Hello world")
    expect_equal(class(valid), "logical")
})

test_that("diamond_is_installed() returns a logical scalar", {
    diamond <- diamond_is_installed()
    expect_equal(class(diamond), "logical")
})

test_that("iqtree_is_installed() returns a logical scalar", {
    iqtree <- iqtree_is_installed()
    expect_equal(class(iqtree), "logical")
})

test_that("species_colors() returns a list for plot_profiles()", {
    species_order <- c(
        "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
        "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
        "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
    ) 
    species_annotation <- data.frame(
        Species = species_order,
        Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
                   "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3), 
                   "Linaceae", "Salicaceae")
    )
    l <- species_colors(species_annotation)
    expect_equal(class(l), "list")
    expect_equal(length(l), 2)
    
})

test_that("valid_seq() checks class of seq object", {
    
    check <- valid_seq(proteomes)
    
    expect_error(
        valid_seq(matrix(NA, 2, 2))
    )
    
    expect_equal(check, TRUE)
})

test_that("valid_annot() checks class of annotation object", {
    
    check <- valid_annot(annotation)
    
    expect_error(
        valid_annot(matrix(NA, 2, 2))
    )
    
    expect_equal(check, TRUE)
})


test_that("valid_blast() checks class of blast list object", {
    
    check <- valid_blast(blast_list)
    
    expect_error(
        valid_blast(matrix(NA, 2, 2))
    )
    
    expect_equal(check, TRUE)
})

