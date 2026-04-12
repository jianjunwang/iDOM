
#' @title Plot Relative Abundance by Group
#' @description Plot relative abundance of molecular data by groups.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data
#' @param mol.group Data frame containing groups for molecular data.
#' @param group_col Column name in mol.group indicating groups. Default: 'group'
#' @return A ggplot object.
#' @rdname plotRA
#' @export 

plotRA <- function(mol.data, mol.group, group_col = "group") {

  # Checking row names consistency
  if(!identical(colnames(mol.data), rownames(mol.group))){
    stop("Mismatch in row names between mol.data and mol.group")
  }
  
  mol.data.RA = t(data.frame(vegan::decostand(mol.data, method = "total"), check.names = FALSE))
  
  merged_data <- cbind(mol.data.RA, group = mol.group[[group_col]])
  
  mean_abundance_by_group <- merged_data %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum)) %>%
    tidyr::pivot_longer(
      cols = -group,
      names_to = "Sample",
      values_to = "Relative_Abundance"
    )

  # Create ggplot object
  p <- ggplot2::ggplot(
    mean_abundance_by_group, 
    ggplot2::aes(x = Sample, y = Relative_Abundance, fill = group)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggplot2::labs(x = "Group", y = "Relative_Abundance") +
    ggplot2::theme_bw () +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  # Return the plot object
  return(p)
}
