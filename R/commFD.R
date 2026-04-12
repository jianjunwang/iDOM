
#' @title Calculate Molecular Functional Diversity
#' @description This function calculates functional diversity indices using molecular data and trait data.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.trait Data frame containing molecular trait data such as mass or other properties associated with each molecule. Default: mol.trait.
#' @param trait.name Character vector giving one or more trait columns in `mol.trait`. Default: "AI_Mod".
#' @param type Character vector specifying which functional diversity indices to return. Available options are `"RaoQ"`, `"FDis"`, `"FEve"`, `"FRic"`, `"FDiv"`, and `"CWM"`. Default: `"RaoQ"`.
#' @return A data frame containing functional diversity indices for each sample.
#' @seealso 
#'  \code{\link[vegan]{decostand}}
#'  \code{\link[fundiversity]{fd_raoq}}
#' @rdname commFD
#' @export
#' @importFrom fundiversity fd_fdis
#' @importFrom fundiversity fd_fdiv
#' @importFrom fundiversity fd_feve
#' @importFrom fundiversity fd_fric
#' @importFrom fundiversity fd_raoq
#' @importFrom vegan decostand

commFD <- function(mol.data, mol.trait, trait.name = c("AI_Mod"), type = c("CWM")){
  
  norm.mol.data = vegan::decostand(mol.data, method = "total")
  type <- unique(as.character(type))
  valid_type <- c("RaoQ", "FDis", "FEve", "FRic", "FDiv", "CWM")
  
  if (length(type) < 1 || any(!type %in% valid_type)) {
    stop("type must contain one or more of: RaoQ, FDis, FEve, FRic, FDiv, CWM.")
  }
  
  if (length(trait.name) < 1 || any(!trait.name %in% colnames(mol.trait))) {
    missing_traits <- setdiff(trait.name, colnames(mol.trait))
    stop(paste("Column(s)", paste(missing_traits, collapse = ", "), "not found in the mol.trait"))
  }
  
  trait_values <- mol.trait[, trait.name, drop = FALSE]
  shared_features <- intersect(colnames(norm.mol.data), rownames(trait_values))
  
  if (length(shared_features) == 0) {
    stop("No shared molecules were found between 'mol.data' columns and 'mol.trait' row names.")
  }
  
  norm.mol.data <- norm.mol.data[, shared_features, drop = FALSE]
  trait_values <- trait_values[shared_features, , drop = FALSE]
  
  if (anyNA(trait_values)) {
    keep <- stats::complete.cases(trait_values)
    norm.mol.data <- norm.mol.data[, keep, drop = FALSE]
    trait_values <- trait_values[keep, , drop = FALSE]
  }
  
  if (nrow(trait_values) == 0) {
    stop("No non-missing trait values are available for the selected trait.")
  }
  
  index_map <- list(
    RaoQ = "fd_raoq",
    FDis = "fd_fdis",
    FEve = "fd_feve",
    FRic = "fd_fric",
    FDiv = "fd_fdiv",
    CWM = "fd_cwm"
  )
  
  result_list <- lapply(type, function(idx) {
    if (idx == "CWM") {
      cwm_values <- as.matrix(norm.mol.data) %*% as.matrix(trait_values)
      cwm_values <- as.data.frame(cwm_values, row.names = rownames(norm.mol.data))
      if (ncol(cwm_values) == 1) {
        colnames(cwm_values) <- "fd_cwm"
      } else {
        colnames(cwm_values) <- paste0("fd_cwm_", colnames(trait_values))
      }
      return(cwm_values)
    }
    
    metric_tbl <- switch(
      idx,
      RaoQ = fundiversity::fd_raoq(traits = trait_values, sp_com = norm.mol.data),
      FDis = fundiversity::fd_fdis(traits = trait_values, sp_com = norm.mol.data),
      FEve = fundiversity::fd_feve(traits = trait_values, sp_com = norm.mol.data),
      FRic = tryCatch(
        fundiversity::fd_fric(traits = trait_values, sp_com = norm.mol.data),
        error = function(e) NULL
      ),
      FDiv = tryCatch(
        fundiversity::fd_fdiv(traits = trait_values, sp_com = norm.mol.data),
        error = function(e) NULL
      )
    )
    
    if (is.null(metric_tbl)) {
      values <- rep(NA_real_, nrow(norm.mol.data))
      names(values) <- rownames(norm.mol.data)
      return(as.numeric(values))
    }
    
    metric_col <- setdiff(colnames(metric_tbl), "site")
    values <- metric_tbl[[metric_col]]
    names(values) <- metric_tbl$site
    values[rownames(norm.mol.data)]
  })
  
  names(result_list) <- type
  
  out <- data.frame(row.names = rownames(norm.mol.data))
  for (idx in type) {
    values <- result_list[[idx]]
    if (is.data.frame(values)) {
      out <- cbind(out, values)
    } else {
      out[[index_map[[idx]]]] <- values
    }
  }
  
  out
}
       