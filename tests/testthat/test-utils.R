
#----Load data------------------------------------------------------------------
data(proteomes)
data(annotation)
data(blast_list)

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
    expect_equal(check, TRUE)
})

test_that("valid_annot() checks class of annotation object", {
    check <- valid_annot(annotation)
    expect_equal(check, TRUE)
})

test_that("valid_blast() checks class of blast list object", {
    check <- valid_blast(blast_list)
    expect_equal(check, TRUE)
})

