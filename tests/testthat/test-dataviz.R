
#----Load data------------------------------------------------------------------
data(network)
data(clusters)
profiles <- phylogenomic_profile(clusters)

#----Start tests----------------------------------------------------------------
test_that("plot_profiles() returns a heatmap", {
    p <- plot_profiles(profiles, cluster_columns = TRUE)
    expect_true("pheatmap" %in% class(p))
})

test_that("plot_network() returns a ggplot object", {
    p1 <- plot_network(network, clusters, cluster_id = 25)
    p2 <- plot_network(network, clusters, cluster_id = clusters$Cluster[1:6])
    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
})
