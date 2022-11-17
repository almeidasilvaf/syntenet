
#' Plot a heatmap of phylogenomic profiles
#'
#' @param profile_matrix A matrix of phylogenomic profiles obtained
#' with \code{phylogenomic_profile}.
#' @param species_annotation A 2-column data frame with species IDs in 
#' the first column (same as column names of profile matrix), and species
#' annotation (e.g., higher-level taxonomic information) in the second column.
#' @param palette A character vector of colors or a character scalar with
#' the name of an RColorBrewer palette. Default: "RdYlBu".
#' @param dist_function Function to use to calculate a distance matrix for
#' synteny clusters. Popular examples include \code{stats::dist},
#' \code{labdsv::dsvdis}, and \code{vegan::vegdist}. Default: stats::dist.
#' @param dist_params A list with parameters to be passed to the function
#' specified in parameter \strong{dist_function}. 
#' Default: list(method = "euclidean").
#' @param clust_function Function to use to cluster the distance matrix
#' returned by the function specified in \code{dist_function}. Examples include
#' \code{stats::hclust} and \code{Rclusterpp::Rclusterpp.hclust}.
#' Default: stats::hclust.
#' @param clust_params A list with additional parameters (if any) to be
#' passed to the function specified in parameter \strong{clust_function}.
#' Default: list(method = "ward.D").
#' @param cluster_species Either a logical scalar (TRUE or FALSE) or 
#' a character vector with the order in which species should be arranged. 
#' TRUE or FALSE indicate whether hierarchical clustering should be applied
#' to rows (species). Ideally, the character vector
#' should contain the order of species in a phylogenetically meaningful way.
#' If users pass a named vector, vector names will be used to rename species.
#' If users have a species tree, they can read it 
#' with \code{treeio::read.tree()}, plot it with \code{ggtree::ggtree()}, 
#' and get the species order from the ggtree object 
#' with \code{ggtree::get_taxa_name()}. Default: FALSE.
#' @param show_colnames Logical indicating whether to show column names (i.e.,
#' cluster IDs) or not. Showing cluster IDs can be useful when visualizing
#' a small subset of them. When visualizing all clusters, cluster IDs are
#' impossible to read. Default: FALSE.
#' @param discretize Logical indicating whether to discretize clusters in
#' 4 categories: 0, 1, 2, and 3+. If FALSE, counts will be
#' log2 transformed. Default: TRUE.
#' @param ... Additional parameters to \code{pheatmap::pheatmap()}.
#'
#' @return A pheatmap object.
#' 
#' @export
#' @rdname plot_profiles
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats dist hclust
#' @examples 
#' data(clusters)
#' profile_matrix <- phylogenomic_profile(clusters)
#' species_order <- c(
#'     "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
#'     "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
#'     "Mnot", "Zjuj", "jcu", "mes", "rco", "lus", "ptr"
#' )
#' species_names <- c(
#'     "V. radiata", "V. angularis", "P. vulgaris", "G. max", "C. cajan",
#'     "T. pratense", "M. truncatula", "A. duranensis", "L. japonicus",
#'     "L. angustifolius", "C. arietinum", "P. mume", "P. persica",
#'     "P. bretschneideri", "M. domestica", "R. occidentalis", 
#'     "F. vesca", "M. notabilis", "Z. jujuba", "J. curcas",
#'     "M. esculenta", "R. communis", "L. usitatissimum", "P. trichocarpa"
#' )
#' names(species_order) <- species_names
#' species_annotation <- data.frame(
#'    Species = species_order,
#'    Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
#'               "Moraceae", "Ramnaceae", rep("Euphorbiaceae", 3), 
#'               "Linaceae", "Salicaceae")
#' )
#' p <- plot_profiles(profile_matrix, species_annotation, 
#'                    cluster_species = species_order)
#'                    
#' p <- plot_profiles(profile_matrix, species_annotation, 
#'                    cluster_species = species_order, 
#'                    discretize = FALSE)
plot_profiles <- function(
        profile_matrix = NULL, 
        species_annotation = NULL,
        palette = "Greens",
        dist_function = stats::dist,
        dist_params = list(method = "euclidean"),
        clust_function = stats::hclust,
        clust_params = list(method = "ward.D"),
        cluster_species = FALSE,
        show_colnames = FALSE,
        discretize = TRUE, ...
) {
    # Set defaults
    breaks <- NA
    legend_labels <- NA
    palette <- colorRampPalette(brewer.pal(7, palette))(100)
    annot_row <- NA
    annot_colors <- NA
    
    # Handle colors and legend for species metadata
    if(is.data.frame(species_annotation)) {
        annot_row <- species_colors(species_annotation)$species_annotation
        annot_colors <- species_colors(species_annotation)$annotation_colors
    }
    
    # Should counts be discretized?
    if(discretize) {
        profile_matrix[profile_matrix >= 3] <- 3
        palette <- c("grey95", "#8FB0D7", "#FDC87A", "#F70C0E")
        breaks <- c(0, 1, 2, 3, 4)
        legend_labels <- c("0", "1", "2", "3+", "")
    } else {
        profile_matrix <- log2(profile_matrix + 1)
    }
    
    # Cluster columns
    ## Get a distance matrix
    dparams <- c(
        list(x = profile_matrix), dist_params
    )
    dist_mat <- do.call(dist_function, dparams)
    ## Cluster based on the distance matrix
    cparams <- c(
        list(d = dist_mat), clust_params
    )
    cluster_columns <- do.call(clust_function, cparams)

    # Handle reordering of species names and change names for better viz
    if(is.character(cluster_species)) {
        if(!all(cluster_species %in% colnames(profile_matrix))) {
            stop("'cluster_species' must match column names of profile matrix.")
        }
        profile_matrix <- profile_matrix[, cluster_species]
        if(!is.null(names(cluster_species))) {
            colnames(profile_matrix) <- names(cluster_species)
            if(is.data.frame(annot_row)) {
                annot_row <- annot_row[cluster_species, , drop = FALSE]
                rownames(annot_row) <- names(cluster_species)
            }
        }
        cluster_species <- FALSE
    }
    # Plot heatmap
    p <- pheatmap::pheatmap(
        t(profile_matrix), color = palette,
        cluster_cols = cluster_columns,
        cluster_rows = cluster_species,
        show_colnames = show_colnames,
        legend_breaks = breaks,
        legend_labels = legend_labels,
        annotation_row = annot_row,
        annotation_colors = annot_colors,
        ...
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
    genes$Group <- as.factor(genes$Cluster) # Attr. 3: Group
    if(is.data.frame(color_by)) {
        genes$Group <- NULL
        names(color_by) <- c("Gene", "Group")
        genes <- merge(genes, color_by, by = "Gene")
        genes$Group <- as.factor(genes$Group)
    }
    graph <- graph_from_data_frame(fedges, directed = FALSE, vertices = genes)
    palette <- custom_palette(1)[seq_len(nlevels(genes$Group))]
    if(nlevels(genes$Group) > 20) {
        palette <- colorRampPalette(custom_palette(1))(nlevels(genes$Group))
    }
    if(interactive) {
        graph_d3 <- networkD3::igraph_to_networkD3(graph, group = genes)
        d <- dim_interactive
        p <- networkD3::forceNetwork(
            Links = graph_d3$links, Nodes = graph_d3$nodes,
            Source = 'source', Target = 'target',
            NodeID = 'name', Group = 'Group', charge = -10,
            height = d[2], width = d[1], Nodesize = 'Degree',
            opacity=0.8, zoom = TRUE, fontSize = 12
        )
    } else {
        n <- ggnetwork::ggnetwork(graph, arrow.gap = 0)
        # Plot graph
        p <- ggplot2::ggplot(n, ggplot2::aes_(x = ~x, y = ~y, 
                                              xend = ~xend, yend = ~yend)) +
            ggnetwork::geom_edges(color = "grey75", alpha = 0.5, 
                                  show.legend = FALSE) +
            ggnetwork::geom_nodes(aes_(size = ~Degree, color = ~Group)) +
            ggplot2::guides(size = "none") +
            ggplot2::scale_color_manual(values = palette) +
            ggnetwork::theme_blank()
    }
    return(p)
}
