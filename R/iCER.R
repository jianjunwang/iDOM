
#' @title The Indicator of Compositional-Level Thermal Responses (iCER)
#' @description The iCER function calculates the Indicator of Compositional-Level Thermal Responses (iCER) based on environmental data (`envi`) and molecular relative abundance data (`data_ra`).
#' @param data_ra A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @param envi A data frame containing environmental variables associated with each sample. The rows should correspond to the samples in `data_ra`, and it should include the variable(s) specified in `vars_list`.
#' @param prop The proportion of samples selected in each permutation for MER (molecule-specific thermal response) calculation. This determines the 'MER dataset' size and affects the stability and variability of the results. Default: 0.8.
#' @param Temp_num The number of temperature levels randomly selected in each permutation.
#' @param end_n The minimum number of times each sample must appear in the 'iCER datasets' during random data partitioning. Increasing this value will result in a greater number of permutations. Default: 10.
#' @param cutoff The minimum proportion of samples in the 'MER dataset' in which a molecule must appear to be retained; molecules below this cutoff are excluded from MER estimation. Default: 0.3.
#' @param vars_list The name of environmental variables (such as 'Temperature') to be used for MER calculation. Default: 'Temperature'.
#' @param significant_only Whether to compute iCER using only molecules with statistically significant MER values (p ≤ `p_value_threshold`). Default: FALSE.
#' @param p_value_threshold The p-value threshold used to assess whether a molecule’s MER is statistically significant. Default: 0.05.
#' @return A list containing the following components:
#' \itemize{
#'   \item `MER`: A data frame of MER values and p-values for molecules in each permutation.
#'   \item `iCER_all_permutation`: A data frame containing iCER, weighted variance (`iCER.var`), standard deviation (`iCER.sd`), and permutation index for each sample, computed from all molecules (if `significant_only = False`).
#'   \item `iCER_sig_permutation`: A data frame containing iCER, weighted variance (`iCER.var`), standard deviation (`iCER.sd`), and permutation index for each sample, computed from molecules with significant MERs (if `significant_only = TRUE`).
#'   \item `iCER`: A data frame containing the aggregated iCER results for each sample across all permutations, including the mean iCER, its weighted variance (`iCER.var`), and standard deviation (`iCER.sd`).
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
#'  \code{\link[Hmisc]{wtd.var}}
#' @rdname iCER
#' @export 
#' @importFrom FD functcomp
#' @importFrom Hmisc wtd.var

iCER <- function(data_ra, envi, prop = 0.80, Temp_num, end_n = 10, cutoff = 0.3, vars_list = "Temperature", significant_only = FALSE, p_value_threshold = 0.05) {
  
  library(dplyr)
  library(tidyverse)
  library(reshape2)
  library(plyr)
  library(vegan)
  library(FD)
  library(Hmisc)
  
  # Prepare data
  envi_tmp <- envi[, vars_list, drop = FALSE]
  
  dat <- tibble(ID = rownames(envi), Temperature = envi[, vars_list])
  
  # Function to generate independent data sets
  random_ID_select <- function(dat, prop, Temp_num, end_n) {
    set.seed(12345)
    num <- round(nrow(dat) * prop)
    Temp_gra_ser <- sort(unique(dat$Temperature))
    
    samp.list <- list()
    unselect.list <- list()
    ID_num <- matrix(nrow = num)
    
    while (TRUE) {
      Temp_gra <- sort(sample(Temp_gra_ser, Temp_num, replace = FALSE))
      random_ID <- sample(which(dat$Temperature %in% Temp_gra), num, replace = FALSE)
      
      # df_check <- dat[random_ID, ]
      # while (!all(Temp_gra %in% df_check$Temperature)) {
      #   random_ID <- sample(which(dat$Temperature %in% Temp_gra), num, replace = FALSE)
      #   df_check <- dat[random_ID, ]
      # }
      
      random_ID <- sort(random_ID)
      
      if (!any(apply(ID_num, 2, function(dup_det) identical(dup_det, random_ID)))) {
        samp.list <- c(samp.list, list(dat[random_ID, 1]))
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
  # peak.rand.list <- c()
  
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
    
    # ## collect peaks
    # peak.go.tmp <- list(peak.go); names(peak.go.tmp) <- kk
    # peak.rand.list <- c(peak.rand.list, peak.go.tmp)
  }
  
  MER.out <- data.frame(MER = rho.out.tmp$spearman.cor, MER.pvalue = rho.out.tmp$spearman.cor.p, Permutation = rho.out.tmp$Permutation, peak = rho.out.tmp$peak)
  
  # (3) Calculate iCER using either ALL molecules with MERs or the molecules with significant MERs
  iCER.out <- data.frame()

  for (k in 1:length(samp.rand.list)) {
    # Filter based on whether significant_only is TRUE or FALSE
    if (significant_only) {
      rho.peak.go0 <- subset(MER.out, Permutation == k & MER.pvalue <= p_value_threshold)   # MERs from "MER dataset" with significant p-values
    } else {
      rho.peak.go0 <- subset(MER.out, Permutation == k)   # MERs from "MER dataset"
    }
    
    rho.peak.go <- rho.peak.go0[, c("MER", "MER.pvalue", "peak")]  
    peak.permu.go <- as.character(rho.peak.go$peak)
    
    unselect.permu.go <- unselect.list[[k]]
    unselect.permu.go <- as.character(unselect.permu.go$ID)
    
    # Relative abundance of molecules in each sample ("iCER dataset")
    data_ra.go <- data_ra[unselect.permu.go, ]  
    data_ra.go <- data_ra.go[, colSums(data_ra.go) > 0]
    
    peak.comm.list <- intersect(peak.permu.go, colnames(data_ra.go))
    
    rho.go <- subset(rho.peak.go, peak %in% peak.comm.list)
    rownames(rho.go) <- rho.go$peak
    rho.go <- rho.go[, "MER", drop=FALSE]
    
    ra.go <- data_ra.go[, rownames(rho.go)]
    ra.go <- ra.go[rowSums(ra.go) > 0, ]
    
    # Compute the iCER for each sample in this permutation
    iCER.tmp <- FD::functcomp(rho.go, as.matrix(ra.go))
    iCER.tmp <- data.frame(ID = rownames(iCER.tmp), iCER = iCER.tmp[, 1])
    
    # Compute the weighted variance and standard deviation of iCER for each sample in this permutation
    if (identical(rownames(rho.go), rownames(t(ra.go)))) {
      x = rho.go[, "MER"]
      wt <- as.data.frame(t(ra.go * 100))
      
      iCER.var.tmp <- sapply(wt, function(wt) wtd.var(x, wt))
      iCER.sd.tmp <- sqrt(iCER.var.tmp)
    }
    
    # Collect iCERs
    iCER.out.tmp <- data.frame(iCER.tmp, iCER.var = iCER.var.tmp, iCER.sd = iCER.sd.tmp, Permutation = k)
    rownames(iCER.out.tmp) = NULL
    iCER.out <- rbind(iCER.out, iCER.out.tmp)
  }
  
  iCER.mean <- ddply(iCER.out, .(ID), summarise, iCER = mean(iCER), iCER.var = mean(iCER.var), iCER.sd = mean(iCER.sd))
  
  if (significant_only) {
    return(list(MER = MER.out, iCER_sig_permutation = iCER.out, iCER = iCER.mean))
  } else {
    return(list(MER = MER.out, iCER_all_permutation = iCER.out, iCER = iCER.mean))
  }
}
