
#' @title Generate Molecular Dendrogram
#' @description This function analyzes molecular traits or transformation data to generate dendrograms. It supports multiple methods for dendrogram generation, including molecular characterization dendrograms (MCD), transformation dendrograms (TD), and transformed weighted characterization dendrograms (TWCD).
#' @param mol.trait A data frame containing molecular traits which must include a column for molecular formula ('MolForm') and may include other specified traits.
#' @param type A character string specifying the type of dendrogram to generate: 'MCD' for Molecular Characterization Dendrogram, 'TD' for Transformation Dendrogram, and 'TWCD' for Transformed Weighted Characterization Dendrogram.
#' @param trait_col A vector of character strings specifying the names of molecular traits to use in the analysis.
#' @param peak.2.peak An optional data frame of peak-to-peak relationships, defaulting to `peak.2.peak`. It is used in the TD and TWCD methods to determine the connections based on transformations. Default: NULL
#' @return A dendrogram object representing the relationships among molecules.
#' @rdname commDendro
#' @export 

commDendro <- function(mol.trait, type, trait_col, peak.2.peak = NULL) {
  
  library(igraph)
  library(phangorn)
  library(vegan)
  library(picante)
  
  options(digits = 10)
  
  # Removing peaks that have no formula assignments
  mol.trait <- mol.trait[!is.na(mol.trait$MolForm),]
  
  # Ensuring that isotopic peaks are removed
  if("C13" %in% colnames(mol.trait)){
    mol.trait <- mol.trait[mol.trait$C13 <= 0,]
  }
  
  # Find column names that are not in mol.trait
  missing_cols <- trait_col[!(trait_col %in% colnames(mol.trait))]
  
  # Check if there are any missing column names
  if (length(missing_cols) > 0) {
    cat("The following columns are not in mol.trait:\n")
    cat(missing_cols, sep = "\n")
  } else {
    cat("All columns are in mol.trait\n")
  }
  
  Mol.Info = mol.trait[, trait_col, drop = F]
  
  if (type == "MCD") {
    
    # Pairwise distance between peaks
    Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info)) # Generating a distance matrix based upon the provided parameters
    
    # Converting the distance matrix into a tree
    tree = as.phylo(hclust(vegdist(Mol.Info, "euclidean"), "average"))
    
  } else if (type == "TD") {
    
    if (is.null(peak.2.peak)) {
      peak.trans <- molTrans(mol.data = mol.data, mol.trait = mol.trait)
      peak.2.peak <- peak.trans$Peak.2.peak %>% 
        select(peak.x, peak.y, Trans.name) %>% 
        rename(from = peak.x, to = peak.y, type = Trans.name) %>% 
        mutate(weight = 1)
    } else {
      peak.2.peak <- peak.2.peak %>% 
        select(peak.x, peak.y, Trans.name) %>% 
        rename(from = peak.x, to = peak.y, type = Trans.name) %>% 
        mutate(weight = 1)
    }

    # Creating the network
    net = graph_from_data_frame(d = peak.2.peak, directed = F)
    
    # The distances command is much better than the similarity measurement
    net.dist = distances(net)
    
    # Finding clusters and determining the distance in the largest
    clus = clusters(net)
    max.clus = which(clus$csize %in% max(clus$csize)) # Finding the largest cluster
    max.clus = names(clus$membership)[which(clus$membership %in% max.clus)] # Finding the members of the largest cluster
 
    net.dist = net.dist[max.clus, max.clus]
    
    # Need to normalize the dissimiarlity to 0-1
    net.dist = (net.dist-min(net.dist))/(max(net.dist)-min(net.dist))
    
    # Generate the UPGMA tree - a neighbor joining one will not work with this large of an object
    tree = as.phylo(hclust(as.dist(net.dist), method = "average"))
    
  } else if (type == "TWCD") {
    
    ### Pairwise molecular distance between peaks
    Mol.Info = as.data.frame(apply(Mol.Info, 2, scale), row.names = row.names(Mol.Info))
    # Generating a distance matrix based upon the provided parameters
    mol.dist = as.matrix(vegdist(Mol.Info, "euclidean"))
    
    if (is.null(peak.2.peak)) {
      peak.trans <- molTrans(mol.data = mol.data, mol.trait = mol.trait)
      peak.2.peak <- peak.trans$Peak.2.peak %>% 
        select(peak.x, peak.y, Trans.name) %>% 
        rename(from = peak.x, to = peak.y, type = Trans.name) %>% 
        mutate(weight = 1)
    } else {
      peak.2.peak <- peak.2.peak %>% 
        select(peak.x, peak.y, Trans.name) %>% 
        rename(from = peak.x, to = peak.y, type = Trans.name) %>% 
        mutate(weight = 1)
    }
    
    ### Determining transformation distance
    # Creating the network
    net = graph_from_data_frame(d = peak.2.peak, directed = F)
    
    # The distances command is much better than the similarity measurement
    net.dist = distances(net)
    
    # Finding clusters and determining the distance in the largest
    clus = clusters(net)
    max.clus = which(clus$csize %in% max(clus$csize)) # Finding the largest cluster
    max.clus = names(clus$membership)[which(clus$membership %in% max.clus)] # Finding the members of the largest cluster
    
    net.dist = net.dist[max.clus, max.clus] # Setting the net dist to that size only; only lost ~4000 peaks by doing this with the dereplicated HJ-Andrews set
    
    # Need to normalize the dissimiarlity to 0-1
    net.dist = (net.dist - min(net.dist))/(max(net.dist) - min(net.dist))
    
    ### Parsing down the data
    q = which(row.names(net.dist) %in% row.names(mol.dist)) # Net dist is generate in the Merged_EdgeNode_Files script - I probably will incorporate it here (or something)
    net.dist.data = net.dist[q, q] # Matching the network distance to the molecular information
    
    q = which(row.names(mol.dist) %in% row.names(net.dist.data))
    mol.dist.data = mol.dist[q, q] # Matching the molecular information to the network distance
    
    # Weighting and tree generation
    weighted.dist = mol.dist.data*net.dist.data # Weighting the tree
    
    # Creating tree
    tree = as.phylo(hclust(as.dist(weighted.dist), method = "average"))
    
  } else {
    stop("Invalid type specified")
  }
  
  return(tree)
}
