
#' @title Plot Van Krevelen Diagram
#' @description Function to plot a Van Krevelen diagram.
#' @param data A data frame containing the H/C ratios and O/C ratios.
#' @param HtoC_ratio The name of the column in `data` containing H/C ratios.
#' @param OtoC_ratio The name of the column in `data` containing O/C ratios.
#' @param title Title of the plot. Default: 'Van Krevelen Diagram'.
#' @return A ggplot object representing the Van Krevelen diagram.
#' @details This function creates a scatter plot of H/C vs O/C ratios.
#' @rdname plotVK
#' @export

plotVK <- function(data, HtoC_ratio, OtoC_ratio, title = "Van Krevelen Diagram") {
  # Check if necessary columns exist
  if (!(HtoC_ratio %in% colnames(data))) {
    stop(paste("Column", HtoC_ratio, "not found in the data"))
  }
  if (!(OtoC_ratio %in% colnames(data))) {
    stop(paste("Column", OtoC_ratio, "not found in the data"))
  }
  
  # Create the base plot
  p <- ggplot(data, aes(x = OtoC_ratio, y = HtoC_ratio))
  
  # Add points and theme
  p <- p + geom_point() +
    labs(x = "O/C Ratio", y = "H/C Ratio") +
    theme_bw()
  
  # Return the plot
  return(p)
}

