
#----Setup----------------------------------------------------------------------
data(network)
data(clusters)

#----Start tests----------------------------------------------------------------
test_that("cluster_network() clusters genes using the infomap algorithm", {
    cluster <- cluster_network(network[1:100, ])
    expect_equal(class(cluster), "data.frame")
    expect_equal(class(cluster$Gene), "character")
    expect_equal(class(cluster$Cluster), "numeric")
    expect_equal(ncol(cluster), 2)
})


test_that("phylogenomic_profile() returns cluster profiles", {
    fclusters <- sample(clusters$Cluster, size = 50, replace = FALSE)
    fclusters <- clusters[clusters$Cluster %in% fclusters, ]
    profiles <- phylogenomic_profile(fclusters)
    expect_true("matrix" %in% class(profiles$profile_matrix))
    expect_equal(ncol(profiles$profile_matrix), 24)
    expect_equal(nrow(profiles$profile_matrix), 
                 length(unique(fclusters$Cluster)))
    expect_equal(class(profiles), "list")
    expect_equal(length(profiles), 2)
    expect_equal(class(profiles$hclust), "hclust")
})
