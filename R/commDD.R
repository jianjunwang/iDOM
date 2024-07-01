
#' @title Calculate Molecular Dendrogram-based Diversity
#' @description This function calculates various Dendrogram-based diversity indices for molecular data.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.dendrogram The dendrogram object representing the relationships among molecules.
#' @param type A character vector specifying the types of diversity indices to calculate. Default: c("DD", "MPD", "MNTD").
#' @return A data frame with the calculated diversity indices for each sample.
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#'  \code{\link[picante]{pd}}, \code{\link[picante]{mpd}}, \code{\link[picante]{mntd}}
#' @rdname commDD
#' @export 
#' @importFrom vegan decostand
#' @importFrom picante pd mpd mntd

commDD <- function(mol.data, mol.dendrogram, type = c("DD", "MPD", "MNTD")){
  library(picante)
  library(vegan)
  
  # Normalize molecular data
  norm.mol.data <- vegan::decostand(mol.data, method = "total")
  
  # Validate the input types
  valid_types <- c("DD", "MPD", "MNTD")
  if (!all(type %in% valid_types)) {
    stop("Invalid type provided. Valid types are: DD, MPD, MNTD.")
  }
  
  # Initialize the result data frame
  comm.DD <- as.data.frame(matrix(data = NA, nrow = nrow(mol.data), ncol = length(type)))
  colnames(comm.DD) <- type
  
  # Define a list of functions for dynamic call
  diversity_functions <- list(
    DD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      return(picante::pd(data.matched, dendrogram.matched)[,"PD"])
    },
    
    MPD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      dendrogram.dist <- cophenetic(dendrogram.matched)
      return(picante::mpd(data.matched, dendrogram.dist, abundance.weighted = T))
    },
    
    MNTD = function(data, dendrogram){
      combined <- match.phylo.comm(dendrogram, data)
      dendrogram.matched <- combined$phy
      data.matched <- combined$comm
      dendrogram.dist <- cophenetic(dendrogram.matched)
      return(picante::mntd(data.matched, dendrogram.dist, abundance.weighted = T)) 
    }
  )
  
  # Iterate over the selected types and calculate the indices
  for (t in type) {
    if (t %in% names(diversity_functions)) {
      comm.DD[[t]] <- diversity_functions[[t]](norm.mol.data, mol.dendrogram)
    }
  }
  
  return(comm.DD)
}
