
#' Proteomes of Ostreococcus sp. species
#'
#' Data obtained from Pico-PLAZA 3.0. Only the translated sequences of primary
#' transcripts were included.
#'
#' @name proteomes
#' @format A list of AAStringSet objects containing 
#' the elements `Olucimarinus`, `Osp_RCC809`, and `Otauri`.
#' @references
#' Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A.,
#' Coppens, F., & Vandepoele, K. (2021). PLAZA 5.0: extending the scope
#' and power of comparative and functional genomics in plants.
#' Nucleic acids research.
#' @examples
#' data(proteomes)
#' @usage data(proteomes)
"proteomes"


#' Genome annotation for Ostreococcus sp. species
#'
#' Data obtained from Pico-PLAZA 3.0. Only annotation data for primary
#' transcripts were included.
#'
#' @name proteomes
#' @format A CompressedGRangesList containing 
#' the elements `Olucimarinus`, `Osp_RCC809`, and `Otauri`.
#' @references
#' Van Bel, M., Silvestri, F., Weitz, E. M., Kreft, L., Botzki, A.,
#' Coppens, F., & Vandepoele, K. (2021). PLAZA 5.0: extending the scope
#' and power of comparative and functional genomics in plants.
#' Nucleic acids research.
#' @examples
#' data(annotation)
#' @usage data(annotation)
"annotation"


#' Synteny network for the legume species G. max, P. vulgaris, and C. cajan
#'
#' @name network
#' @format An edgelist (i.e., a 2-column data frame with node 1 in 
#' column 1 and node 2 in column 2).
#' @references
#' Zhao, T., & Schranz, M. E. (2019). Network-based microsynteny analysis 
#' identifies major differences and genomic outliers in mammalian and 
#' angiosperm genomes. Proceedings of the National Academy of 
#' Sciences, 116(6), 2165-2174.
#' @examples
#' data(network)
#' @usage data(network)
"network"



#' Synteny network clusters for the legume species 
#' G. max, P. vulgaris, and C. cajan
#'
#' @name clusters
#' @format A 2-column data frame containing the following variables:
#' \describe{
#'   \item{Gene}{Gene ID}
#'   \item{Cluster}{Cluster ID}
#' }
#' @references
#' Zhao, T., & Schranz, M. E. (2019). Network-based microsynteny analysis 
#' identifies major differences and genomic outliers in mammalian and 
#' angiosperm genomes. Proceedings of the National Academy of 
#' Sciences, 116(6), 2165-2174.
#' @examples
#' data(clusters)
#' @usage data(clusters)
"clusters"
