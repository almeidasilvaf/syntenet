
#' Plot a heatmap of phylogenomic profiles
#'
#' @param profiles A list of phylogenomic profiles as returned
#' by \code{phylogenomic_profile()}.
#' @param species_annotation A 2-column data frame with species IDs in 
#' the first column (same as column names of profile matrix), and species
#' annotation (e.g., higher-level taxonomic information) in the second column.
#' @param palette A character vector of colors or a character scalar with
#' the name of an RColorBrewer palette. Default: "RdYlBu".
#' @param cluster_species Either a logical scalar (TRUE or FALSE) or 
#' a character vector with the order of columns. Ideally, the character vector
#' should contain the order of species in a phylogenetically meaningful way.
#' If users have a species tree, they can read it 
#' with \code{treeio::read.tree()} and get the species order 
#' with \code{ggtree::get_taxa_name()}.
#' @param cluster_columns Either a logical scalar (TRUE or FALSE) or an hclust
#' object. Default: hclust object in the \strong{profiles} list as returned
#' by \code{phylogenomic_profile()}.
#' @param discretize Logical indicating whether to discretize clusters in
#' 4 categories: 0, 1, 2, and 3+. 
#'
#' @return A pheatmap object.
#' 
#' @export
#' @rdname plot_profiles
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette 
#' @importFrom RColorBrewer brewer.pal
#' @examples 
#' data(clusters)
#' profiles <- phylogenomic_profile(clusters)
#' species_order <- c(
#'     "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
#'     "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
#'     "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
#' ) 
#' species_annotation <- data.frame(
#'    Species = species_order,
#'    Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
#'               "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3), 
#'               "Linaceae", "Salicaceae")
#')
#' p <- plot_profiles(profiles, species_annotation, 
#'                    cluster_species = species_order, 
#'                    cluster_columns = TRUE)
plot_profiles <- function(
        profiles = NULL, 
        species_annotation = NULL,
        palette = "Greens",
        cluster_species = FALSE,
        cluster_columns = profiles$hclust,
        discretize = TRUE
) {
    # Set defaults
    breaks <- NA
    legend_labels <- NA
    profile_matrix <- t(profiles$profile_matrix)
    palette <- colorRampPalette(brewer.pal(7, palette))(100)
    annot_row <- NA
    annot_colors <- NA
    if(is.data.frame(species_annotation)) {
        annot_row <- species_colors(species_annotation)$species_annotation
        annot_colors <- species_colors(species_annotation)$annotation_colors
    }
    
    if(discretize) {
        profile_matrix[profile_matrix >= 3] <- 3
        palette <- c("grey95", "#8FB0D7", "#FDC87A", "#F70C0E")
        breaks <- c(0, 1, 2, 3, 4)
        legend_labels <- c("0", "1", "2", "3+", "")
    } else {
        profile_matrix <- log2(profile_matrix + 1)
    }
    if(is.character(cluster_species)) {
        if(!all(cluster_species %in% rownames(profile_matrix))) {
            stop("'cluster_species' must match column names of profile matrix.")
        }
        profile_matrix <- profile_matrix[cluster_species, ]
        cluster_species <- FALSE
    }
    
    p <- pheatmap::pheatmap(
        profile_matrix, color = palette,
        cluster_cols = cluster_columns,
        cluster_rows = cluster_species,
        show_colnames = FALSE,
        legend_breaks = breaks,
        legend_labels = legend_labels,
        annotation_row = annot_row,
        annotation_colors = annot_colors
    )
    return(p)
}


#' Plot network
#' 
#' @param network The synteny network represented as an edge list, which is
#' a 2-column data frame with each member of the anchor pair in a column.
#' @param clusters A 2-column data frame with the variables \strong{Gene} 
#' and \strong{Cluster} representing gene ID and cluster ID, respectively,
#' exactly as returned by \code{cluster_network}.
#' @param cluster_id Character scalar or vector with cluster ID. If more than one 
#' cluster is passed as input, clusters are colored differently. 
#' @param color_by Either "cluster" or a 2-column data frame with 
#' gene IDs in the first column and variable to be used for 
#' coloring (e.g., taxonomic information) in the second column.
#' @param interactive Logical scalar indicating whether to display an 
#' interactive network or not. Default: FALSE.
#' @param dim_interactive Numeric vector of length 2 with the window 
#' dimensions of the interactive plot. If \strong{interactive} is set to FALSE,
#' this parameter is ignored.
#' 
#' @return A ggplot object with the network.
#' 
#' @importFrom ggnetwork ggnetwork theme_blank geom_edges geom_nodes
#' @importFrom networkD3 igraph_to_networkD3 forceNetwork
#' @importFrom igraph graph_from_data_frame
#' @importFrom ggplot2 aes_ ggplot scale_color_manual guides
#' @importFrom grDevices colorRampPalette
#' @import intergraph
#' @export
#' @rdname plot_network
#' @examples 
#' data(network)
#' data(clusters)
#' # Option 1: 1 cluster
#' cluster_id <- 25
#' plot_network(network, clusters, cluster_id)
#' 
#' # Option 2: 2 clusters
#' cluster_id <- c(25, 1089)
#' plot_network(network, clusters, cluster_id)
#' 
#' # Option 3: custom annotation for coloring
#' species_order <- c(
#'     "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
#'     "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
#'     "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
#' ) 
#' species_annotation <- data.frame(
#'    Species = species_order,
#'    Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
#'               "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3), 
#'               "Linaceae", "Salicaceae")
#')
#' genes <- unique(c(network$node1, network$node2))
#' gene_df <- data.frame(
#'     Gene = genes,
#'     Species = unlist(lapply(strsplit(genes, "_"), head, 1))
#' )
#' gene_df <- merge(gene_df, species_annotation)[, c("Gene", "Family")]
#' 
#' plot_network(network, clusters, cluster_id = 25, color_by = gene_df)
plot_network <- function(network = NULL, clusters = NULL, 
                         cluster_id = NULL, color_by = "cluster",
                         interactive = FALSE, 
                         dim_interactive = c(600, 600)) {
    
    requireNamespace("intergraph", quietly=TRUE)
    names(network) <- c("node1", "node2")
    #----Handle node attributes------------------------------------------------
    genes <- clusters[clusters$Cluster %in% cluster_id, ] # Att. 1: Cluster
    fedges <- network[network$node1 %in% genes$Gene, ]
    fedges <- fedges[fedges$node2 %in% genes$Gene, ]
    deg <- as.data.frame(table(c(fedges$node1, fedges$node2)))
    names(deg) <- c("Gene", "Degree")
    genes <- merge(genes, deg) # Attr. 2: Degree
    genes$Class <- as.factor(genes$Cluster) # Attr. 3: Class
    if(is.data.frame(color_by)) {
        genes$Class <- as.factor(color_by[, 2][color_by[, 1] %in% genes$Gene])
    }
    graph <- graph_from_data_frame(fedges, directed = FALSE, vertices = genes)
    palette <- custom_palette(1)[seq_len(nlevels(genes$Class))]
    if(nlevels(genes$Class) > 20) {
        palette <- colorRampPalette(custom_palette(1))(nlevels(genes$Class))
    }
    if(interactive) {
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = genes)
        d <- dim_interactive
        p <- networkD3::forceNetwork(
            Links = graph_d3$links, Nodes = graph_d3$nodes,
            Source = 'source', Target = 'target',
            NodeID = 'name', Group = 'Class',
            height = d[2], width = d[1], Nodesize = 'Degree',
            opacity=0.8, zoom = TRUE, fontSize = 12
        )
    } else {
        n <- ggnetwork::ggnetwork(graph, arrow.gap = 0)
        # Plot graph
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, xend = ~xend, yend = ~yend)) +
            ggnetwork::geom_edges(color = "grey75", alpha = 0.5, show.legend=FALSE) +
            ggnetwork::geom_nodes(ggplot2::aes_(size = ~Degree, color = ~Class)) +
            ggplot2::guides(size = "none") +
            ggplot2::scale_color_manual(values = palette) +
            ggnetwork::theme_blank()
    }
    return(p)
}




