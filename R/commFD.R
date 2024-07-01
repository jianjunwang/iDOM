
#' @title Calculate Molecular Functional Diversity
#' @description This function calculates functional diversity indices using molecular data and trait data.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.trait Data frame containing molecular trait data such as mass or other properties associated with each molecule. Default: mol.trait.
#' @param trait.name The name of the trait column in the `mol.trait`. Default: "AI_Mod".
#' @return A data frame containing functional diversity indices for each sample.
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#'  \code{\link[fundiversity]{fd_raoq}}
#' @rdname commFD
#' @export 
#' @importFrom vegan decostand
#' @importFrom fundiversity fd_raoq

commFD <- function(mol.data, mol.trait, trait_col = "AI_Mod"){
  library(fundiversity)
  library(vegan)
  
  # Normalize mol.data using total method
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  
  # Check if necessary columns exist
  if (!(trait_col %in% colnames(mol.trait))) {
    stop(paste("Column", HtoC_ratio, "not found in the mol.trait"))
  }
  
  # Calculate functional diversity using the specified trait column
  comm.FD = fundiversity::fd_raoq(sp_com = norm.mol.data, traits = mol.trait[, trait_col, drop = F])
  
  return(comm.FD)
}
