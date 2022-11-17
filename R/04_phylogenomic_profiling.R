

#' Cluster the synteny network using the Infomap algorithm
#'
#' @param network A network represented as an edge list, which is a 2-column
#' data frame with node 1 in the first column and node 2 in the second column.
#' In a synteny network, node 1 and node are the anchor pairs.
#' @param clust_function Function to be used to cluster the network. It must
#' be one the functions from the \strong{cluster_*} family in 
#' the \strong{igraph} package (e.g., cluster_infomap, cluster_leiden, etc). 
#' Default: igraph::cluster_infomap.
#' @param clust_params A list with additional parameters (if any) to be passed 
#' to the igraph clustering function. Default: NULL (no additional parameters).
#' 
#' @return A 2-column data frame with the following variables:
#' \describe{
#'   \item{Gene}{Gene ID.}
#'   \item{Cluster}{Cluster ID as identified by infomap.}
#' }
#'
#' @importFrom igraph graph_from_data_frame cluster_infomap simplify
#' @export
#' @rdname cluster_network
#' @examples 
#' data(network)
#' clusters <- cluster_network(network[1:500, ])
cluster_network <- function(network = NULL, 
                            clust_function = igraph::cluster_infomap,
                            clust_params = NULL) {
    
    # Create a graph and remove loops
    graph <- igraph::graph_from_data_frame(network, directed = FALSE)
    graph <- igraph::simplify(graph)

    # Cluster the graph    
    cparams <- c(
        list(graph = graph), clust_params
    )
    clusters <- do.call(clust_function, cparams)
    
    # Create a data frame of genes and their corresponding clusters
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
#' @return A matrix of i rows and j columns containing 
#' the number of genes in cluster i for each species j. The number of rows 
#' is equal to the number of clusters in \strong{clusters}, and 
#' the number of columns is equal to the number of species 
#' in \strong{clusters}.
#'
#' @export
#' @rdname phylogenomic_profile
#' @examples 
#' data(clusters)
#' profiles <- phylogenomic_profile(clusters)
phylogenomic_profile <- function(clusters = NULL) {
    
    # Add species info
    clusters$Species <- vapply(
        strsplit(clusters$Gene, "_"), `[`, 1, FUN.VALUE = character(1)
    )
    # Create a profile matrix
    profile_matrix <- table(clusters$Cluster,clusters$Species)
    profile_matrix <- matrix(
        profile_matrix, ncol = ncol(profile_matrix), 
        dimnames = dimnames(profile_matrix)
    )
    
    return(profile_matrix)
}


#' Find group-specific clusters based on user-defined species classification
#'
#' @param profile_matrix A matrix of phylogenomic profiles obtained
#' with \code{phylogenomic_profile}.
#' @param species_annotation A 2-column data frame with species IDs in 
#' the first column (same as column names of profile matrix), and species
#' annotation (e.g., higher-level taxonomic information) in the second column.
#' @param min_percentage Numeric scalar with the minimum percentage of species
#' in a group to consider group specificity. For instance,
#' if a given cluster is present in only 1 group of species, but in less
#' than \strong{min_percentage} of the species for this group, it will not
#' be considered a group-specific cluster. This filtering criterion is useful 
#' to differentiate group-specific clusters (e.g., family-specific) from 
#' subgroup-specific clusters (e.g., genus-specific). Default: 50.
#' 
#' @return A data frame with the following variables:
#' \describe{
#'   \item{Group}{To which group of species the cluster is specific.}
#'   \item{Percentage}{Percentage of species from the group that are
#'   represented by the cluster.}
#'   \item{Cluster}{Cluster ID.}
#' }
#' @importFrom stats reshape
#' @export
#' @rdname find_GS_clusters
#' @examples
#' data(clusters)
#' profile_matrix <- phylogenomic_profile(clusters)
#' 
#' # Species annotation
#' species_order <- c(
#'     "vra", "van", "pvu", "gma", "cca", "tpr", "mtr", "adu", "lja",
#'     "Lang", "car", "pmu", "ppe", "pbr", "mdo", "roc", "fve",
#'     "Mnot", "Zjuj", "hlu", "jcu", "mes", "rco", "lus", "ptr"
#' )
#' species_annotation <- data.frame(
#'    Species = species_order,
#'    Family = c(rep("Fabaceae", 11), rep("Rosaceae", 6),
#'               "Moraceae", "Ramnaceae", "Cannabaceae",
#'                rep("Euphorbiaceae", 3), "Linaceae", "Salicaceae")
#' )
#' gs_clusters <- find_GS_clusters(profile_matrix, species_annotation)
find_GS_clusters <- function(profile_matrix = NULL, 
                             species_annotation = NULL,
                             min_percentage = 50) {
    if(!is.data.frame(species_annotation)) {
        stop("Species annotation must be a 2-column data frame.")
    }
    diff <- setdiff(colnames(profile_matrix), species_annotation$Species)
    if(length(diff) != 0) {
        diffname <- paste0(diff, collapse = "\n")
        message("Could not find annotation for species: \n", diffname)
    }
    names(species_annotation) <- c("Species", "Group")
    freq_by_group <- as.data.frame(table(species_annotation$Group))
    n <- length(unique(species_annotation$Group))
    
    # Reshape to long format
    bmat <- as.data.frame(profile_matrix)
    bmat[bmat > 0] <- 1
    bmat_long <- reshape(bmat, idvar = "Cluster", ids = rownames(bmat),
                         times = names(bmat), timevar = "Species", 
                         varying = list(names(bmat)), direction = "long")
    names(bmat_long) <- c("Species", "Presence", "Cluster")
    
    # Add group information
    prof_grouped <- merge(bmat_long, species_annotation)
    
    # Look for group-specific clusters
    prof_list <- split(prof_grouped, prof_grouped$Cluster)
    gs_df <- Reduce(rbind, lapply(seq_along(prof_list), function(x) {
        gs <- NULL
        filt_df <- prof_list[[x]][prof_list[[x]]$Presence == 1, ]
        group_count <- as.data.frame(table(filt_df$Group))
        group_count <- merge(freq_by_group, group_count, by = "Var1", 
                             all.x = TRUE)
        group_count$perc <- (group_count$Freq.y / group_count$Freq.x) * 100
        group_count[is.na(group_count)] <- 0
        
        # How many groups have 0% of species?
        perc_zero <- sum(group_count$perc == 0)
        if(perc_zero == n - 1) {
            # How many families are represented by >x% of the genes?
            fgroup_count <- group_count[group_count$perc > min_percentage, ]
            if(nrow(fgroup_count) == 1) {
                gs <- fgroup_count[, c("Var1", "perc")]
                gs$Cluster <- names(prof_list)[x]
            }
        }
        return(gs)
    }))
    
    if(is.data.frame(gs_df)) {
        names(gs_df) <- c("Group", "Percentage", "Cluster")
        gs_df$Percentage <- round(gs_df$Percentage, 2)
    }
    return(gs_df)
}
