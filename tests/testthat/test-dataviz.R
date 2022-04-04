
#----Load data------------------------------------------------------------------
data(clusters)
profiles <- phylogenomic_profile(clusters)

#----Start tests----------------------------------------------------------------
test_that("plot_profiles() returns a heatmap", {
    p <- plot_profiles(profiles, cluster_columns = TRUE)
    expect_true("pheatmap" %in% class(p))
})
