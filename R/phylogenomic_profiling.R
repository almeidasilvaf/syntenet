

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
#' library(cogeqc)
#' data(synnet)
#' network <- synnet[1:500, ]
#' clusters <- cluster_network(network)
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

