
#' @title Molecular Grouping Analysis
#' @description This function groups molecular data based on specific traits.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data
#' @param mol.trait Data frame containing molecular trait data such as mass or other properties associated with each molecule. Default: mol.trait
#' @param group The grouping criteria, Default: 'Reactivity_Activity'
#' @param Trans_threshold The threshold for transformations, Default: c(1, 10)
#' @return A data frame with grouped molecules.
#' @rdname molGroup
#' @export 

molGroup <- function(mol.data, mol.trait, group = "Reactivity_Activity", HtoC_ratio, Peak.trans.profile = NULL, Trans_threshold = c(1,10)) {
  
  library(dplyr)
  
  # Check if the column names of mol.data match the row names of mol.trait
  if(!identical(colnames(mol.data), rownames(mol.trait))){
    stop("Mismatch in names: The column names of mol.data do not match the row names of mol.trait. Please check the consistency of your data.")
  }
  
  if (group == "Reactivity_Activity") {
    # Reactivity: labile: H/C ≥ 1.5; recalcitrant: H/C < 1.5
    
    # Check if necessary columns exist
    if (!(HtoC_ratio %in% colnames(mol.trait))) {
      stop(paste("Column", HtoC_ratio, "not found in the mol.trait"))
    }
    
    mol.HC = mol.trait %>% 
      select(HtoC_ratio) %>% 
      mutate(peak = rownames(.)) %>% 
      relocate(peak, .before = HtoC_ratio)
    
    if(is.null(Peak.trans.profile)){
      Transformation.results <- molTrans(mol.data = mol.data, mol.trait = mol.trait)
      peak.profile <- Transformation.results[[2]]
    } else {
      peak.profile <- Peak.trans.profile
    }
    
    dim(peak.profile)
    
    mol.trait.group <- peak.profile %>% 
      select(peak, num.trans.involved.in) %>% 
      left_join(mol.HC, by = "peak") %>% 
      rename(HtoC_ratio = colnames(.)[3]) %>% 
      mutate(group = case_when(HtoC_ratio >= 1.5 & num.trans.involved.in > Trans_threshold[2] ~"Labile_Active",
                               HtoC_ratio >= 1.5 & num.trans.involved.in <= Trans_threshold[1] ~"Labile_Inactive",
                               HtoC_ratio < 1.5 & num.trans.involved.in > Trans_threshold[2] ~"Recalcitrant_Active",
                               HtoC_ratio < 1.5 & num.trans.involved.in <= Trans_threshold[1] ~"Recalcitrant_Inactive"))
    
  }else {
    # Check if necessary columns exist
    if (!(group %in% colnames(mol.trait))) {
      stop(paste("Column", group, "not found in the mol.trait"))
    }
    
    mol.trait.group <- mol.trait %>% 
      select(group) %>% 
      mutate(peak = rownames(.)) %>% 
      relocate(peak, .before = group)
  }
  
  return(mol.trait.group)
}
