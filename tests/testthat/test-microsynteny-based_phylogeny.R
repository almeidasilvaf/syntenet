
#----Load data------------------------------------------------------------------
data(clusters)
profile_matrix <- phylogenomic_profile(clusters)$profile_matrix
set.seed(123)

#----Start tests----------------------------------------------------------------
test_that("binarize_and_transpose() binarizes and transposes the profiles", {
    tmat <- binarize_and_transpose(profile_matrix)
    expect_true("matrix" %in% class(tmat))
    expect_equal(nrow(profile_matrix), ncol(tmat))
    expect_equal(ncol(profile_matrix), nrow(tmat))
    expect_equal(length(which(tmat > 1)), 0)
})

test_that("profiles2phylip() returns path to a PHYLIP file", {
    tmat <- binarize_and_transpose(profile_matrix)
    path <- profiles2phylip(tmat)
    expect_equal(class(path), "character")
    expect_equal(length(path), 1)
})

test_that("infer_microsynteny_phylogeny() infers a phylogeny", {
    tmat <- binarize_and_transpose(profile_matrix)
    # Leave only some legumes and P. mume as an outgroup
    included <- c("gma", "pvu", "vra", "van", "cca", "pmu")
    tmat <- tmat[rownames(tmat) %in% included, ]
    
    # Remove non-variable sites
    tmat <- tmat[, colSums(tmat) != length(included)]
    
    phylo <- character(length = 10)
    if(iqtree_is_installed()) {
        phylo <- infer_microsynteny_phylogeny(tmat, outgroup = "pmu", 
                                              threads = 1)
    }
    expect_equal(class(phylo), "character")
    expect_equal(length(phylo), 10)
})
