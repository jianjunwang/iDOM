
#' @title The Indicator of Molecular Dark Matter Effects (iDME)
#' @description The function calculates the iDME by measuring the percentage change in the average value of a specified network metric between the “KK” (known-known) and “DK” (dark-known) networks.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.dark.matter Data frame containing molecular dark matter where rows represent different samples or conditions, and columns represent individual unclassified molecules. Default: mol.dark.matter.
#' @param occu.rate The threshold for retaining DOM molecules observed in a specified proportion of the total samples. Default is 0.3.
#' @param bootstrap The number of 'KK' and 'DK' networks to be randomly generated for each sample; defaults to 100 (more might be better, less probably not).
#' @param Network.size The size of 'KK' and 'DK' networks. Default: 400.
#' @param sparcc.R An optional precomputed SparCC correlation matrix. If NULL, the function will compute it. Default is NULL.
#' @param sparcc.R.threshold A numeric value specifying the threshold for the absolute value of the SparCC correlation matrix to filter out uncorrelated or weakly correlated interactions. Correlations with an absolute value less than 0.3 (i.e., -0.3 < correlation < 0.3) will be removed. Default is 0.3.
#' @return A data frame with calculated iDME, iDME.intra, and iDME.inter for each sample.
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#'  \code{\link[SpiecEasi]{sparcc}}
#'  \code{\link[igraph]{graph_from_adjacency_matrix}}, \code{\link[igraph]{V}}, \code{\link[igraph]{degree}}
#'  \code{\link[bipartite]{specieslevel}}
#'  \code{\link[dplyr]{rename}}
#' @rdname iDME
#' @export 
#' @importFrom vegan decostand
#' @importFrom SpiecEasi sparcc
#' @importFrom igraph graph_from_adjacency_matrix V degree
#' @importFrom bipartite specieslevel
#' @importFrom dplyr rename

iDME <- function(mol.data, mol.dark.matter, occu.rate = 0.3, bootstrap = 100, Network.size = 400, sparcc.R = NULL, sparcc.R.threshold = 0.3){

  # Load necessary libraries
  library(vegan)
  library(SpiecEasi)
  library(dplyr)
  library(bipartite)
  library(igraph)
  
  # Checking row names consistency
  if(!identical(rownames(mol.data), rownames(mol.dark.matter))){
    stop("Row names of mol.data and mol.dark.matter do not match.")
  }else{
    mol.data.all <- cbind(mol.data, mol.dark.matter)
    mol.data.all.RA <- vegan::decostand(mol.data.all, method = "total")
    mol.data.all.round <- as.matrix(round(100000 * mol.data.all.RA, 0))
  }

  mol.data.all.occu <- mol.data.all.round[, colSums(mol.data.all.round > 0) >= nrow(mol.data.all.round) * occu.rate];dim(mol.data.all.occu)

  if (is.null(sparcc.R)) {
    
    sparcc.R <- SpiecEasi::sparcc(data = mol.data.all.occu)
    sparcc.R <- sparcc.R$Cor
    colnames(sparcc.R) <- rownames(sparcc.R) <- colnames(mol.data.all.occu)
    
    sparcc.R.thresh <- ifelse(abs(sparcc.R) >= sparcc.R.threshold, abs(sparcc.R), 0)
    
  }else{
    
    # Check consistency of column and row names
    if (!identical(colnames(mol.data.all.occu), rownames(sparcc.R))) {
      stop("Column names of mol.data.all.occu do not match row names of sparcc.R")
    }
    if (!identical(colnames(mol.data.all.occu), colnames(sparcc.R))) {
      stop("Column names of mol.data.all.occu do not match column names of sparcc.R")
    }
    
    sparcc.R.thresh <- ifelse(abs(sparcc.R) >= sparcc.R.threshold, abs(sparcc.R), 0)
  }
  
  Net_Degree <- data.frame()
  
  ii = 1; jj = 1
  for (ii in 1:nrow(mol.data.all.occu)) {
    mol.data.sample <- mol.data.all.occu[ii, , drop = F];dim(mol.data.sample)
    mol.data.sample <- mol.data.sample[, colSums(mol.data.sample) > 0, drop = F]; dim(mol.data.sample)

    mol.known.occu <- mol.data.sample[, colnames(mol.data.sample) %in% colnames(mol.data), drop = F];dim(mol.known.occu)
    mol.unknown.occu <- mol.data.sample[, colnames(mol.data.sample) %in% colnames(mol.dark.matter), drop = F];dim(mol.unknown.occu)

    for (jj in 1:bootstrap) {
      
      # Network nodes
      Known.seq <- sample(ncol(mol.known.occu), Network.size, replace = FALSE)
      K1.seq <- sample(Known.seq, Network.size/2, replace = FALSE)
      K2.seq <- Known.seq[!Known.seq %in% K1.seq]

      Unknown.seq <- sample(ncol(mol.unknown.occu), Network.size/2, replace = FALSE)

      Net.K1 = mol.known.occu[, K1.seq, drop = F]; dim(Net.K1)
      Net.K2 = mol.known.occu[, K2.seq, drop = F]; dim(Net.K2)
      Net.D = mol.unknown.occu[, Unknown.seq, drop = F]; dim(Net.D)

      # Known networks
      Known_Net.all = sparcc.R.thresh[c(colnames(Net.K1), colnames(Net.K2)), c(colnames(Net.K1), colnames(Net.K2))];dim(Known_Net.all)
      Known.all <- igraph::graph_from_adjacency_matrix(as.matrix(Known_Net.all), mode = "undirected", weighted = TRUE, diag = FALSE)
      Known.all.degree <- data.frame(
        nodes_id = igraph::V(Known.all)$name,Degree = igraph::degree(Known.all)
      );rownames(Known.all.degree) = NULL

      Net_Known.Inner.K1 <- Known_Net.all[colnames(Net.K1), colnames(Net.K1)];dim(Net_Known.Inner.K1)
      Net_Known.Inner.K2 <- Known_Net.all[colnames(Net.K2), colnames(Net.K2)];dim(Net_Known.Inner.K2)

      Known.Inner.K1 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_Known.Inner.K1), mode = "undirected", weighted = TRUE, diag = FALSE)
      Known.Inner.K2 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_Known.Inner.K2), mode = "undirected", weighted = TRUE, diag = FALSE)

      Known.Inner.degree <- data.frame(
        nodes_id = c(igraph::V(Known.Inner.K1)$name,igraph::V(Known.Inner.K2)$name),Degree = c(igraph::degree(Known.Inner.K1),igraph::degree(Known.Inner.K2))
      );rownames(Known.Inner.degree) = NULL

      Net_Known.Inter <- Known_Net.all[colnames(Net.K2), colnames(Net.K1)];dim(Net_Known.Inter)
      Known.Inter.higher <- bipartite::specieslevel(Net_Known.Inter, index="degree", level = "higher");dim(Known.Inter.higher)
      Known.Inter.lower <- bipartite::specieslevel(Net_Known.Inter, index="degree", level = "lower");dim(Known.Inter.lower)
      Known.Inter.degree = data.frame(nodes_id = c(rownames(Known.Inter.higher),rownames(Known.Inter.lower)), Degree = c(Known.Inter.higher$degree,Known.Inter.lower$degree))

      if (identical(Known.Inter.degree$nodes_id,Known.all.degree$nodes_id) == F) {
        Known.Inter.degree.0 = data.frame(nodes_id = Known.all.degree$nodes_id[!Known.all.degree$nodes_id %in% Known.Inter.degree$nodes_id],Degree = 0)
        Known.Inter.Degree = rbind(Known.Inter.degree,Known.Inter.degree.0)
      }else{
        Known.Inter.Degree = Known.Inter.degree
      }

      Net_Known_degree = Known.all.degree %>%
        full_join(Known.Inner.degree,by = "nodes_id") %>%
        full_join(Known.Inter.Degree,by = "nodes_id") %>%
        dplyr::rename(KK_all.degree = Degree.x, KK_inner.degree = Degree.y, KK_inter.degree = Degree) %>%
        select(KK_all.degree, KK_inner.degree, KK_inter.degree) %>%
        summarise_all(mean, na.rm = TRUE)

      # DK_Networks -------------------------------------------------------------
      DK_Net.all = sparcc.R.thresh[c(colnames(Net.K1), colnames(Net.D)), c(colnames(Net.K1), colnames(Net.D))];dim(DK_Net.all)

      DK.all <- igraph::graph_from_adjacency_matrix(as.matrix(DK_Net.all), mode = "undirected", weighted = TRUE, diag = FALSE)
      DK.all.degree <- data.frame(
        nodes_id = igraph::V(DK.all)$name,Degree = igraph::degree(DK.all)
      );rownames(DK.all.degree) = NULL

      Net_DK.Inner.K1 <- DK_Net.all[colnames(Net.K1), colnames(Net.K1)];dim(Net_DK.Inner.K1)
      Net_DK.Inner.D <- DK_Net.all[colnames(Net.D), colnames(Net.D)];dim(Net_DK.Inner.D)

      DK.Inner.K1 <- igraph::graph_from_adjacency_matrix(as.matrix(Net_DK.Inner.K1), mode = "undirected", weighted = TRUE, diag = FALSE)
      DK.Inner.D <- igraph::graph_from_adjacency_matrix(as.matrix(Net_DK.Inner.D), mode = "undirected", weighted = TRUE, diag = FALSE)

      DK.Inner.degree <- data.frame(
        nodes_id = c(igraph::V(DK.Inner.K1)$name,igraph::V(DK.Inner.D)$name),Degree = c(igraph::degree(DK.Inner.K1),igraph::degree(DK.Inner.D))
      );rownames(DK.Inner.degree) = NULL

      # Inter
      Net_DK.Inter <- DK_Net.all[colnames(Net.D), colnames(Net.K1)];dim(Net_DK.Inter)
      DK.Inter.higher <- bipartite::specieslevel(Net_DK.Inter, index="degree", level = "higher");dim(DK.Inter.higher)
      DK.Inter.lower <- bipartite::specieslevel(Net_DK.Inter, index="degree", level = "lower");dim(DK.Inter.lower)
      DK.Inter.degree = data.frame(nodes_id = c(rownames(DK.Inter.higher),rownames(DK.Inter.lower)), Degree = c(DK.Inter.higher$degree,DK.Inter.lower$degree))

      if (!identical(DK.Inter.degree$nodes_id, DK.all.degree$nodes_id)) {
        DK.Inter.degree.0 = data.frame(nodes_id = DK.all.degree$nodes_id[!DK.all.degree$nodes_id %in% DK.Inter.degree$nodes_id], Degree = 0)
        DK.Inter.Degree = rbind(DK.Inter.degree, DK.Inter.degree.0)
      }else{
        DK.Inter.Degree = DK.Inter.degree
      }

      Net_DK_degree = DK.all.degree %>%
        full_join(DK.Inner.degree,by = "nodes_id") %>%
        full_join(DK.Inter.Degree,by = "nodes_id") %>%
        dplyr::rename(DK_all.degree = Degree.x, DK_inner.degree = Degree.y, DK_inter.degree = Degree) %>%
        select(DK_all.degree, DK_inner.degree, DK_inter.degree) %>%
        summarise_all(mean, na.rm = TRUE)

      Net_degree = cbind(Net_Known_degree, Net_DK_degree, bootstrap = jj, sample = rownames(mol.data.sample))

      if (ii == 1 & jj == 1) {
        Net_Degree = Net_degree
      }else{
        Net_Degree = rbind(Net_Degree, Net_degree)
      }
      print(ii)
    }
  }
  
  iDME.sample = Net_Degree %>% 
    group_by(sample) %>% 
    summarise(across(starts_with("KK_"), mean, na.rm = TRUE),
              across(starts_with("DK_"), mean, na.rm = TRUE)) %>% 
    mutate(iDME = (DK_all.degree/KK_all.degree) - 1,
           iDME.intra = (DK_inner.degree-KK_inner.degree)/KK_all.degree,
           iDME.inter = (DK_inter.degree-KK_inter.degree)/KK_all.degree) %>% 
    select(sample, iDME, iDME.intra, iDME.inter)
    
    return(iDME.sample)
  
}
