

#' Cluster the synteny network using the Infomap algorithm
#'
#' @param network A network represented as an edge list, which is a 2-column
#' data frame with node 1 in the first column and node 2 in the second column.
#' In a synteny network, node 1 and node are the anchor pairs.
#' 
#' 
#' @return A 2-column data frame with the following variables:
#' \describe{
#'   \item{Gene}{Gene ID.}
#'   \item{Cluster}{Cluster ID as identified by infomap.}
#' }
#'
#'
#' @importFrom igraph graph_from_data_frame cluster_infomap simplify
#' @export
#' @rdname cluster_network
#' @examples 
#' data(network)
#' clusters <- cluster_network(network[1:500, ])
cluster_network <- function(network = NULL) {
    
    graph <- igraph::graph_from_data_frame(network, directed = FALSE)
    graph <- igraph::simplify(graph)
    
    clusters <- igraph::cluster_infomap(graph)
    clusters_df <- data.frame(
        Gene = clusters$names,
        Cluster = clusters$membership
    )
    return(clusters_df)
}


#' Perform phylogenomic profiling for synteny network clusters
#'
#' @param clusters A 2-column data frame with variables \strong{Gene} and
#' \strong{Cluster} as returned by \code{cluster_network}.
#' 
#' @return A matrix of i rows and j columns containing the number of genes 
#' in cluster i for each species j. The number of rows is equal to 
#' the number of clusters in \strong{clusters}, and the number of columns 
#' is equal to the number of species in \strong{clusters}.
#'
#' @importFrom vegan vegdist 
#' @importFrom stats hclust
#' @export
#' @rdname phylogenomic_profile
#' @examples 
#' data(clusters)
#' profiles <- phylogenomic_profile(clusters)
phylogenomic_profile <- function(clusters = NULL) {
    
    # Add species info and create profile matrix
    clusters$Species <- vapply(strsplit(clusters$Gene, "_"), `[`, 1, 
                               FUN.VALUE = character(1))
    profile_matrix <- table(clusters$Cluster,clusters$Species)
    profile_matrix <- matrix(profile_matrix, ncol = ncol(profile_matrix), 
                             dimnames = dimnames(profile_matrix))
    
    # Calculate matrix of Jaccard distances
    dist_mat <- vegan::vegdist(log2(profile_matrix + 1), method = "jaccard")
    
    # Cluster with ward.D based on Jaccard distances
    clust_mat <- stats::hclust(dist_mat, method = "ward.D")
    
    # Reorder rows based on clustering
    fprofile_matrix <- profile_matrix[clust_mat$order, ]
    return(fprofile_matrix)
}

