
#' @title The specialization index H2' of DOM-microbe associations
#' @description This function calculates the network-level specialization of all interacting trophic levels in DOM-microbe bipartite networks, including full, negative and positive networks.
#' @param mol.data Sample/DOM matrix with samples in the rows and DOM molecules in the columns (compositional/abundant data).
#' @param micro.data Sample/Microbe matrix with samples in the rows and Microbial species in the columns (compositional/abundant data).
#' @param occu.rate The threshold for retaining bacterial species or DOM molecules observed in a specified proportion of the total samples. Default is 0.3.
#' @param sparcc.R An optional precomputed SparCC correlation matrix. If NULL, the function will compute it. Default: NULL
#' @param sparcc.R.threshold The threshold of including correlations between DOM molecules and bacterial species in bipartite networks. Default: 0.3
#' @param N Number of null models to be generated; defaults to 100 (more might be better, less probably not).
#' @param Null.model Null model type. Can be given as an integer or name: 1/"r2dtable", 2/"swap.web", 3/"vaznull", 4/"shuffle.web"; allows for partial match of names. Default: 'shuffle.web'.
#' @return Returns a data.frame, which contains standardised H2', Observed H2', P-value and Network type. Standardised: Standardised H2'. It is standardised by using a null modelling approach. 
#' Observed: Observed H2'. It ranges between 0 (complete generalization) and 1 (complete specialization). p.value: The significance of differences between observed and random H2'. Network.type: Network type. Can be given as a name: Full (full network), 
#' Negative (negative network), Positive (positive network).
#' @details H2' is a network-level property to describe how much the two trophic levels are interacting with each other in a bipartite network. 
#' For example, H2' is used to quantify the specialization of DOM-microbe associations at a network level. 
#' Specifically, elevated H2' values convey that there is a high degree of specialization between DOM and microbes. 
#' By contrast, lower H2' values reflect a more generalized bipartite network where different DOM molecules can be used by a large range of bacterial taxa.
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#'  \code{\link[SpiecEasi]{sparcc}}
#'  \code{\link[bipartite]{nullmodel}}, \code{\link[bipartite]{networklevel}}
#' @rdname H2
#' @export 
#' @importFrom vegan decostand
#' @importFrom SpiecEasi sparcc
#' @importFrom bipartite nullmodel networklevel

H2 <- function(mol.data, micro.data, occu.rate = 0.3, sparcc.R = NULL, sparcc.R.threshold = 0.3, N = 100, Null.model = "shuffle.web") {
  
  library(vegan)
  library(bipartite)
  library(SpiecEasi)
  
  # Standardize and round microbial and DOM data
  norm.mol.data <- vegan::decostand(mol.data, method = "total")
  norm.mol.data.round <- as.matrix(round(100000 * norm.mol.data, 0))
  
  norm.micro.data <- vegan::decostand(micro.data, method = "total")
  norm.micro.data.round <- as.matrix(round(100000 * norm.micro.data, 0))

  # Checking row names consistency
  if(!identical(rownames(norm.mol.data.round), rownames(norm.micro.data.round))){
    stop("Row names of mol.data and micro.data do not match.")
  }else{
    mol.micro.data <- cbind(norm.mol.data.round, norm.micro.data.round)
  }
  
  # Filter columns based on occurrence threshold
  mol.micro.data.keep <- mol.micro.data[, colSums(mol.micro.data > 0) >= nrow(mol.micro.data) * occu.rate]
  dim(mol.micro.data.keep)
  
  if (is.null(sparcc.R)) {
    
    sparcc.R <- SpiecEasi::sparcc(data = mol.micro.data.keep)
    sparcc.R <- sparcc.R$Cor
    colnames(sparcc.R) <- rownames(sparcc.R) <- colnames(mol.micro.data.keep)
    
    sparcc.R.thresh <- ifelse(abs(sparcc.R) >= sparcc.R.threshold, sparcc.R, 0)
    
  }else{
    
    # Check consistency of column and row names
    if (!identical(colnames(mol.micro.data.keep), rownames(sparcc.R))) {
      stop("Column names of mol.data.all.occu do not match row names of sparcc.R")
    }
    if (!identical(colnames(mol.micro.data.keep), colnames(sparcc.R))) {
      stop("Column names of mol.data.all.occu do not match column names of sparcc.R")
    }
    
    sparcc.R <- as.matrix(sparcc.R)
    sparcc.R.thresh <- ifelse(abs(sparcc.R) >= sparcc.R.threshold, sparcc.R, 0)
  }
  
  # Extract relevant subset of adjusted correlation matrix
  sparcc.R.thresh_DOM_Micro <- sparcc.R.thresh[rownames(sparcc.R.thresh) %in% colnames(mol.data), colnames(sparcc.R.thresh) %in% colnames(micro.data)]
  DOM_Micro.int <- round(100000 * sparcc.R.thresh_DOM_Micro, 0)
  
  # Perform network analysis
  classes <- c("Full", "Positive", "Negative")
  Index <- NULL
  
  for (class in classes) {
    if (class == "Full") {
      adj.cor.int = ifelse(abs(DOM_Micro.int) > 0, abs(DOM_Micro.int), 0)
    }else if(class == "Positive"){
      adj.cor.int = ifelse(DOM_Micro.int > 0, DOM_Micro.int, 0) 
    }else{
      adj.cor.int = ifelse(DOM_Micro.int < 0, abs(DOM_Micro.int), 0) 
    }
    
    adj.cor.int.class <- adj.cor.int[rowSums(adj.cor.int) > 0, colSums(adj.cor.int) > 0];dim(adj.cor.int.class)
    
    nulls <- bipartite::nullmodel(adj.cor.int.class, N = N, method = Null.model)

    weighted.indices = c("H2")
    
    for (jj in 1:nrow(mol.micro.data.keep)) {
      site_data <- mol.micro.data.keep[jj, , drop = FALSE];dim(site_data)
      site_data <- site_data[,colSums(site_data) > 0,drop = F];dim(site_data)
      
      site.Micro <- site_data[, colnames(site_data) %in% colnames(micro.data), drop = F];dim(site.Micro)
      site.DOM <- site_data[, colnames(site_data) %in% colnames(mol.data), drop = F];dim(site.DOM)
      
      adj.cor.int.site = adj.cor.int.class[rownames(adj.cor.int.class) %in% colnames(site.DOM),colnames(adj.cor.int.class) %in% colnames(site.Bac)];dim(adj.cor.int.site)
      
      nulls.site = list()
      
      for (i in 1:N) {
        colnames(nulls[[i]]) = colnames(adj.cor.int.class)
        rownames(nulls[[i]]) = rownames(adj.cor.int.class)
        
        nulls.site[[i]] <- nulls[[i]][rownames(adj.cor.int.site),colnames(adj.cor.int.site)];dim(nulls.site[[i]])
      }
      
      index.obs = bipartite::networklevel(adj.cor.int.site, weighted = T, index = weighted.indices)
      index.null = unlist(sapply(nulls.site, bipartite::networklevel, weighted = T, index = weighted.indices))
      
      index.rand.mean = mean(index.null)
      index.rand.sd = sd(index.null)
      index.ses = (index.obs - index.rand.mean) / index.rand.sd
      praw = sum(index.null > index.obs) / length(index.null)
      index.obs.p = ifelse(praw > 0.5, 1-praw, praw)
      
      index <- data.frame(Index = names(index.obs), Observed = index.obs, Standardised = index.ses, p.value = index.obs.p, Network.type = class, Sample = rownames(site_data))
      rownames(index) = NULL
      Index <- rbind(Index, index)
      
      if (jj %% 5 == 0) {
        message(paste("Sample:", rownames(site_data), "; Null.model.type:", class))
      }
      
    }
  }
  
  return(Index)
}





