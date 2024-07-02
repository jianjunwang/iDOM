
#' @title Quantifying Assembly Processes of DOM assemblages
#' @description This function quantifies the relative influences of deterministic and stochastic processes governing the assembly of DOM assemblages.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.dendrogram The dendrogram object representing the relationships among molecules.
#' @param mol.group Optional data frame containing group information for molecular data. Default: NULL.
#' @param group_col Column name in `mol.group` indicating groups. This parameter is required if `mol.group` is provided.
#' @return A list with the following elements: 
#' \item{ratio}{A data frame showing the overall percentage of turnovers governed by each process.}
#' \item{result}{A matrix showing betaMNTD, Bray-Curtis dissimilarity, betaNTI, RC, and the identified process governing each turnover. Each row represents a turnover between two samples.}
#' @rdname commProc
#' @export 
#' @importFrom picante match.phylo.comm
#' @importFrom dplyr filter
#' @importFrom rlang sym

commProc <- function(mol.data, mol.dendrogram, mol.group = NULL, group_col = 'group'){
  
  library(iCAMP)
  library(picante)
  library(dplyr)
  library(rlang)
  
  # Helper function to process a group
  assembly_process <- function(data, dendrogram) {
    dd <- cophenetic(dendrogram)
    qp <- qpen(comm = data, pd = dd, rand.time = 1000, ab.weight = FALSE)
    list(ratio = data.frame(ratio = qp$ratio), result = data.frame(result = qp$result))
  }
  
  results <- list()
  
  if (is.null(mol.group)) {
    phylo <- picante::match.phylo.comm(phy = mol.dendrogram, comm = mol.data)
    res <- assembly_process(phylo$comm, phylo$phy)
    
    results <- list(ratio = res$ratio, result = res$result)
    
  } else {
    groups <- unique(mol.group[[group_col]])
    
    for (group in groups) {
      mol.group.filtered <- mol.group %>% dplyr::filter(!!rlang::sym(group_col) == !!group)
      mol.data.group <- mol.data[, colnames(mol.data) %in% mol.group.filtered$peak]
      
      phylo <- picante::match.phylo.comm(phy = mol.dendrogram, comm = mol.data.group)
      res <- assembly_process(phylo$comm, phylo$phy)
      
      results[[length(results) + 1]] <- list(ratio = res$ratio, result = res$result)
      names(results)[length(results)] <- group
    }
  }
  return(results)
}