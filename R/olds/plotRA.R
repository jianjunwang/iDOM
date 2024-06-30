
#' @title Plot Relative Abundance by Group
#' @description Plot relative abundance of molecular data by groups.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data
#' @param mol.group Data frame containing groups for molecular data.
#' @param group_col Column name in mol.group indicating groups. Default: 'group'
#' @return A ggplot object.
#' @rdname plotRA
#' @export 

plotRA <- function(mol.data, mol.group, group_col = "group") {
  library(ggplot2)
  library(dplyr)
  library(vegan)
  
  # Checking row names consistency
  if(!identical(colnames(mol.data), rownames(mol.group))){
    stop("Mismatch in row names between mol.data and mol.group")
  }
  
  mol.data.RA = t(data.frame(decostand(mol.data, method = "total"), check.names = F))
  
  merged_data <- cbind(mol.data.RA, mol.group[, group_col, drop = F])
  
  mean_abundance_by_group <- merged_data %>%
    group_by(group) %>%
    summarise(across(where(is.numeric), sum)) %>% 
    gather(2:ncol(.), key = "Sample", value = "Relative_Abundance")

  # Create ggplot object
  p <- ggplot(mean_abundance_by_group, aes(x = Sample, y = Relative_Abundance, fill = group)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Group", y = "Relative_Abundance") +
    theme_bw () +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Return the plot object
  return(p)
}
