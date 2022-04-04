
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