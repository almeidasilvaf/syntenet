
#----Setup----------------------------------------------------------------------
data(network)

#----Start tests----------------------------------------------------------------
test_that("cluster_network() clusters genes using the infomap algorithm", {
    cluster <- cluster_network(network[1:100, ])
    expect_equal(class(cluster), "data.frame")
    expect_equal(class(cluster$Gene), "character")
    expect_equal(class(cluster$Cluster), "numeric")
    expect_equal(ncol(cluster), 2)
})
