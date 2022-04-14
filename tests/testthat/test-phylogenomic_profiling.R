
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

test_that("find_GS_clusters() finds group-specific clusters", {
    profile_matrix <- phylogenomic_profile(clusters)$profile_matrix
    
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
    gs_clusters <- find_GS_clusters(profile_matrix, species_annotation)
    
    expect_equal(class(gs_clusters), "data.frame")
    expect_equal(names(gs_clusters), c("Group", "Percentage", "Cluster"))
})

