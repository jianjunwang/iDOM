
#' @title Calculate Molecular Taxonomic Diversity
#' @description This function calculates various taxonomic diversity indices for molecular data.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param type A character vector specifying the diversity indices to calculate. Options are "Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", and "Chao1". Default: Richness.
#' @return A data frame containing calculated diversity indices for each sample.
#' @seealso 
#'  \code{\link[vegan]{decostand}}, \code{\link[vegan]{diversity}}, \code{\link[vegan]{specpool}}
#' @rdname commTD
#' @export 
#' @importFrom vegan decostand specnumber diversity estimateR

commTD <- function(mol.data, type = c("Richness")) {
  library(vegan)
  
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  
  # Validate the input types
  valid_types <- c("Richness", "Shannon", "Simpson", "Invsimpson", "PielouEven", "Chao1")
  if (!all(type %in% valid_types)) {
    stop("Invalid type provided. Valid types are: Richness, Shannon, Simpson, Invsimpson, PielouEven, and Chao1.")
  }
  
  # Initialize the result data frame
  comm.TD <- as.data.frame(matrix(data = NA, nrow = nrow(mol.data), ncol = length(type)))
  colnames(comm.TD) <- type
  
  # Define a list of functions for dynamic call
  diversity_functions <- list(
    Richness = function(data) vegan::specnumber(data),
    Shannon = function(data) vegan::diversity(data, index = "shannon"),
    Simpson = function(data) vegan::diversity(data, index = "simpson"),
    Invsimpson = function(data) vegan::diversity(data, index = "invsimpson"),
    PielouEven = function(data) vegan::diversity(data, index = "shannon") / log(vegan::specnumber(data)),
    Chao1 = function(data) {
      norm_mol_counts_100k = round(norm.mol.data * 100000)
      return(vegan::estimateR(norm_mol_counts_100k, method = "chao")[2,])
    }
  )
  
  # Iterate over the selected types and calculate the indices
  for (t in type) {
    if (t %in% names(diversity_functions)) {
      comm.TD[[t]] <- diversity_functions[[t]](norm.mol.data)
    }
  }
  
  return(comm.TD)
}
