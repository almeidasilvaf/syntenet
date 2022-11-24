
#----Setup----------------------------------------------------------------------
data(network)
data(clusters)
set.seed(123)

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
    
    expect_true("matrix" %in% class(profiles))
    expect_equal(nrow(profiles), length(unique(fclusters$Cluster)))
})

test_that("find_GS_clusters() finds group-specific clusters", {
    
    profile_matrix <- phylogenomic_profile(clusters)
    
    species_order <- c(
        "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
        "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
        "Mnot", "Zjuj", "hlu", "jcu", "mes", "rco", "lus", "ptr"
    )
    species_annotation <- data.frame(
       Species = species_order,
       Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
                  "Moraceae", "Ramnaceae", "Cannabaceae",
                   rep("Euphorbiaceae", 3), "Linaceae", "Salicaceae")
    )
    
    gs_clusters <- find_GS_clusters(profile_matrix, species_annotation)
    
    expect_error(
        find_GS_clusters(profile_matrix, matrix(NA, 2, 2))
    )
    
    expect_equal(class(gs_clusters), "data.frame")
    expect_equal(names(gs_clusters), c("Group", "Percentage", "Cluster"))
})

