
#' @title The Indicator of Compositional-Level Thermal Responses (iCER)
#' @description The iCER function calculates the Indicator of Compositional-Level Thermal Responses (iCER) based on environmental data (`envi`) and molecular relative abundance data (`data_ra`).
#' @param data_ra A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @param envi A data frame containing environmental variables associated with each sample. The rows should correspond to the samples in `data_ra`, and it should include the variable(s) specified in `vars_list`.
#' @param prop The proportion of the dataset to be used for MER (Molecule-specific thermal responses) calculation. This controls how much of the original dataset is sampled for analysis. Default: 0.8.
#' @param Temp_num The number of temperature levels randomly selected in each permutation.
#' @param end_n Ensures that each sample appears at least this many times in the iCER datasets during random data partitioning. Increasing this value will result in a greater number of permutations. Default: 10.
#' @param cutoff The threshold for retaining molecules observed in a specified proportion of the samples within the MER dataset. Molecules present in less than this proportion of samples will be excluded from the MER calculation. Default: 0.3.
#' @param vars_list The name or list of environmental variables (such as temperature) to be used in the correlation analysis. Default: 'Temperature'.
#' @param significant_only A boolean indicating whether to calculate iCER using only molecules with significant MERs. If set to TRUE, only molecules with significant correlations (based on `p_value_threshold`) are included. Default: FALSE.
#' @param p_value_threshold The significance level for determining whether a molecule’s correlation is considered significant. Default: 0.05
#' @return A list containing the following components:
#' \itemize{
#'   \item `MER`: A data frame of MER values for molecules.
#'   \item `CER_all_out`: A data frame of iCER values calculated from all molecules (if `significant_only = FALSE`).
#'   \item `CER_sig_out`: A data frame of iCER values calculated from molecules with significant MERs (if `significant_only = TRUE`).
#'   \item `CER_mean`: A data frame of the mean iCER value for each sample.
#' }
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  # Example of calculating iCER using all molecules
#'  result_all <- iCER(data_ra, envi, significant_only = FALSE)
#'  
#'  # Example of calculating iCER using only significant molecules
#'  result_sig <- iCER(data_ra, envi, significant_only = TRUE, p_value_threshold = 0.05)
#'  }
#' }
#' @seealso 
#'  \code{\link[FD]{functcomp}}
#' @rdname iCER
#' @export 
#' @importFrom FD functcomp

iCER <- function(data_ra, envi, prop = 0.80, Temp_num, end_n = 10, cutoff = 0.3, vars_list = "Temperature", significant_only = FALSE, p_value_threshold = 0.05) {
  
  library(dplyr)
  library(tidyverse)
  library(reshape2)
  library(plyr)
  library(vegan)
  library(FD)
  
  # Prepare data
  envi$cDOM.IDs <- rownames(envi)
  envi_tmp <- envi[, vars_list, drop = FALSE]
  
  dat <- tibble(ID = envi$cDOM.IDs, Temperature = envi[, vars_list])
  
  # Function to generate independent data sets
  random_ID_select <- function(dat, prop, Temp_num, end_n) {
    set.seed(12345)
    num <- round(nrow(dat) * prop)
    ele_gra_ser <- sort(unique(dat$Temperature))
    
    samp.list <- list()
    unselect.list <- list()
    ID_num <- matrix(nrow = num)
    
    while (TRUE) {
      ele_gra <- sort(sample(ele_gra_ser, Temp_num, replace = FALSE))
      random_ID <- sample(which(dat$Temperature %in% ele_gra), num, replace = FALSE)
      
      # df_check <- dat[random_ID, ]
      # while (!all(ele_gra %in% df_check$Temperature)) {
      #   random_ID <- sample(which(dat$Temperature %in% ele_gra), num, replace = FALSE)
      #   df_check <- dat[random_ID, ]
      # }
      
      random_ID <- sort(random_ID)
      
      if (!any(apply(ID_num, 2, function(dup_det) identical(dup_det, random_ID)))) {
        samp.list <- c(samp.list, list(dat[random_ID, ]))
        unselect.list <- c(unselect.list, list(dat[-random_ID, 1]))
      }
      
      ID_num <- cbind(ID_num, random_ID)
      df <- do.call(cbind, unselect.list)
      
      if (all(sapply(dat$ID, function(x) sum(x == df)) >= end_n)) {
        break 
      }
    }
    ID_num <- ID_num[,-1]
    
    return(list(samp.list, unselect.list))
  }
  
  # (1) Generate independent data sets ------------------------------------------
  rand.id.tmp <- random_ID_select(dat, prop, Temp_num, end_n)
  samp.rand.list <- rand.id.tmp[[1]]
  unselect.list <- rand.id.tmp[[2]]
  
  # (2) Calculate MER using "MER dataset" ---------------------------------------
  rho.out.tmp <- data.frame()
  peak.rand.list <- c()
  
  for (kk in 1:length(samp.rand.list)) {
    samp.list.go <- samp.rand.list[[kk]]
    samp.list.go <- as.character(samp.list.go$ID)
    
    data_ra_go2 <- data_ra[samp.list.go, ]
    data_ra_go2 <- data_ra_go2[, colSums(data_ra_go2) > 0]
    data_pa_go2 <- decostand(data_ra_go2, method = "pa")
    data_ra_go2_keep <- data_ra_go2[, colSums(data_pa_go2) >= nrow(data_pa_go2) * cutoff]
    
    ## combining temperature
    dat.tmp <- merge(envi_tmp, data_ra_go2_keep, by.x="row.names", by.y="row.names", all.y = TRUE)
    dat.tmp.long <- reshape2::melt(dat.tmp, id = c(colnames(dat.tmp)[c(1:(1+length(vars_list)))]))
    
    ## peak
    peak.go <- colnames(data_ra_go2_keep)
    
    for (ii in 1:length(peak.go)) {
      dat.tmp.long.i <- subset(dat.tmp.long, variable %in% peak.go[ii])
      
      if (nrow(dat.tmp.long.i) < 5) next
      if (nrow(dat.tmp.long.i) >= 5) {
        tmp <- dat.tmp.long.i[, c("Temperature", "value")]
        corr.out.tmp <- data.frame(spearman.cor = with(tmp, cor(value, Temperature, method="spearman")), 
                                   spearman.cor.p = with(tmp, cor.test(value, Temperature, method = "spearman"))$p.value,
                                   Permutation = kk,
                                   peak = peak.go[ii]
        )
        
        ## collect MERs
        rho.out.tmp <- rbind(rho.out.tmp, corr.out.tmp)
      }
    }
    
    ## collect peaks
    peak.go.tmp <- list(peak.go); names(peak.go.tmp) <- kk
    peak.rand.list <- c(peak.rand.list, peak.go.tmp)
  }
  
  MER.out <- data.frame(MER = rho.out.tmp$spearman.cor, MER.pvalue = rho.out.tmp$spearman.cor.p, Permutation = rho.out.tmp$Permutation, peak = rho.out.tmp$peak)
  
  # (3) Calculate iCER using either ALL molecules with MERs or the molecules with significant MERs
  CER.out <- data.frame()
  
  for (k in 1:length(peak.rand.list)) {   
    # Filter based on whether significant_only is TRUE or FALSE
    if (significant_only) {
      rho.peak.go0 <- subset(rho.out.tmp, Permutation == k & spearman.cor.p <= p_value_threshold)   # MERs from "MER dataset" with significant p-values
    } else {
      rho.peak.go0 <- subset(rho.out.tmp, Permutation == k)   # MERs from "MER dataset"
    }
    
    rho.peak.go <- rho.peak.go0[, c("spearman.cor", "spearman.cor.p", "peak")]  
    peak.permu.go <- as.character(rho.peak.go$peak)
    
    unselect.permu.go <- unselect.list[[k]]
    unselect.permu.go <- as.character(unselect.permu.go$ID)
    
    data_ra.go <- data_ra[unselect.permu.go, ]  # Relative abundance of molecules in each sample ("iCER dataset")
    data_ra.go <- data_ra.go[, colSums(data_ra.go) > 0]
    
    peak.comm.list <- intersect(peak.permu.go, colnames(data_ra.go))
    
    rho.go <- subset(rho.peak.go, peak %in% peak.comm.list)  
    rownames(rho.go) <- rho.go$peak
    rho.go <- rho.go[, "spearman.cor", drop=FALSE]
    
    ra.go <- data_ra.go[, rownames(rho.go)]   
    ra.go <- ra.go[rowSums(ra.go) > 0, ]
    
    CER.tmp <- FD::functcomp(rho.go, as.matrix(ra.go))   # iCER calculation
    CER.tmp2 <- data.frame(CER.tmp, cDOM.IDs = rownames(CER.tmp))
    
    # Collect iCERs
    CER.out.tmp <- data.frame(CER.tmp2, Permutation = k)
    CER.out <- rbind(CER.out, CER.out.tmp)
  }
  
  CER.mean <- ddply(CER.out, .(cDOM.IDs), summarise, iCER=mean(spearman.cor))
  
  colnames(CER.out)[1:2] <- c("iCER","ID")
  rownames(CER.out) <- NULL
  colnames(CER.mean)[1] <- "ID"
  
  if (significant_only) {
    return(list(MER = MER.out, CER_sig_out = CER.out, CER_mean = CER.mean))
  } else {
    return(list(MER = MER.out, CER_all_out = CER.out, CER_mean = CER.mean))
  }
}
