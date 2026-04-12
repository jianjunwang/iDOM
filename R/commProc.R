
#' @title Quantifying Assembly Processes of DOM assemblages
#' @description This function quantifies the relative influences of deterministic and stochastic processes governing the assembly of DOM assemblages.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.dendrogram The dendrogram object representing the relationships among molecules.
#' @param mol.group Optional data frame containing group information for molecular data. Default: NULL.
#' @param group_col Column name in `mol.group` indicating groups. This parameter is required if `mol.group` is provided.
#' The row names of `mol.group` should identify the molecule columns in `mol.data`.
#' @return A list with the following elements: 
#' \item{ratio}{A data frame showing the overall percentage of turnovers governed by each process.}
#' \item{result}{A matrix showing betaMNTD, Bray-Curtis dissimilarity, betaNTI, RC, and the identified process governing each turnover. Each row represents a turnover between two samples.}
#' @rdname commProc
#' @references A Hu, K-S Jang, F Meng, J Stegen, A J Tanentzap, M Choi, J T Lennon, J Soininen, and J Wang. 2022. 
#' Microbial and Environmental Processes Shape the Link between Organic Matter Functional Traits and Composition. 
#' *Environmental Science & Technology*. [https://doi.org/10.1021/acs.est.2c01432](https://pubs.acs.org/doi/abs/10.1021/acs.est.2c01432); [PDF](https://jjwang.name/data/uploads/publications/Hu-et-al-2022-EST.pdf)
#' 
#' F Meng, A Hu, K-S Jang, and J Wang. 2025. 
#' iDOM: Statistical analysis of dissolved organic matter characterized by high-resolution mass spectrometry. 
#' *mLife*. [https://onlinelibrary.wiley.com/doi/10.1002/mlf2.70002](https://onlinelibrary.wiley.com/doi/10.1002/mlf2.70002)
#' 
#' @export 
#' @importFrom picante match.phylo.comm
#' @importFrom dplyr filter
#' @importFrom rlang sym

commProc <- function(mol.data, mol.dendrogram, mol.group = NULL, group_col = 'group'){

  # Helper function to process a group
  assembly_process <- function(data, dendrogram) {
    dd <- stats::cophenetic(dendrogram)
    qp <- iCAMP::qpen(comm = data, pd = dd, rand.time = 1000, ab.weight = FALSE)
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
      mol.group.filtered <- mol.group %>% 
        dplyr::filter(!!rlang::sym(group_col) == !!group)
      mol.data.group <- mol.data[, colnames(mol.data) %in% rownames(mol.group.filtered), drop = FALSE]
      
      phylo <- picante::match.phylo.comm(phy = mol.dendrogram, comm = mol.data.group)
      res <- assembly_process(phylo$comm, phylo$phy)
      
      results[[length(results) + 1]] <- list(ratio = res$ratio, result = res$result)
      names(results)[length(results)] <- group
    }
  }
  return(results)
}
