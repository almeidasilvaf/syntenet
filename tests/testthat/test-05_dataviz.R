
#----Load data------------------------------------------------------------------
data(network)
data(clusters)
profiles <- phylogenomic_profile(clusters)

## Data frame of species annotation
species_order <- c(
    "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
    "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
    "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
)
species_names <- c(
    "V. radiata", "V. angularis", "P. vulgaris", "G. max", "C. cajan",
    "T. pratense", "M. truncatula", "A. duranensis", "L. japonicus",
    "L. angustifolius", "C. arietinum", "P. mume", "P. persica",
    "P. bretschneideri", "M. domestica", "R. occidentalis",
    "F. vesca", "M. notabilis", "Z. jujuba", "J. curcas",
    "M. esculenta", "R. communis", "L. usitatissimum", "P. trichocarpa"
)
names(species_order) <- species_names
species_annotation <- data.frame(
   Species = species_order,
   Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
              "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3),
              "Linaceae", "Salicaceae")
)

# color_by data frame
genes <- unique(c(network$node1, network$node2))
gene_df <- data.frame(
    Gene = genes,
    Species = unlist(lapply(strsplit(genes, "_"), head, 1))
)
gene_df <- merge(gene_df, species_annotation)[, c("Gene", "Family")]
gene_df$Family[1:20] <- paste0("Random", 1:20)
    
#----Start tests----------------------------------------------------------------
test_that("plot_profiles() returns a heatmap", {
    
    p1 <- plot_profiles(profiles, discretize = FALSE)
    p2 <- plot_profiles(
        profiles, 
        species_annotation = species_annotation,
        cluster_species = species_order,
        discretize = TRUE
    )
    p3 <- plot_profiles(
        profiles, 
        species_annotation <- species_annotation[-1, ],
        discretize = FALSE
    )
    
    cspecies <- c(species_order[1:5], "fake_species")
    expect_error(
        plot_profiles(
            profiles, 
            species_annotation = species_annotation,
            cluster_species = cspecies,
            discretize = TRUE
        )
    )
    
    expect_true("pheatmap" %in% class(p1))
    expect_true("pheatmap" %in% class(p2))
    expect_true("pheatmap" %in% class(p3))
})


test_that("plot_network() returns a ggplot object", {
    
    p1 <- plot_network(
        network, 
        clusters, 
        cluster_id = c(
            227, 470, 277, 1059, 759, 685, 102, 210, 196, 48, 895, 635, 
            1097, 521, 1297, 303, 873, 258, 1175, 1212
        ), 
        color_by = gene_df
    )
    p2 <- plot_network(network, clusters, cluster_id = clusters$Cluster[1:6])
    p3 <- plot_network(network, clusters, cluster_id = 1, interactive = TRUE)
    
    expect_true("ggplot" %in% class(p1))
    expect_true("ggplot" %in% class(p2))
    expect_true("forceNetwork" %in% class(p3))
})
