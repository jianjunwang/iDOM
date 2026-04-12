
#' @title The Indicator of Molecular Habitat Affinity (iHaf)
#' @description The iHaf function calculates the Indicator of Molecular Habitat Affinity (iHaf) based on molecular relative abundance data of habitat1 (data.ra.habi1) and habitat2 (data.ra.habi2).
#' @param data.ra.habi1 A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @param data.ra.habi2 A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @return A list containing the following components:
#' \describe{
#'   \item{iHaf}{A data frame of iHaf values for molecules.}
#'   \item{iHaf_wm_habi1}{A data frame of iHaf_wm values for samples of habitat1.}
#'   \item{iHaf_wm_habi2}{A data frame of iHaf_wm values for samples of habitat2.}
#' }
#' @references Y Cui, A Hu, J Stegen, and J Wang. 2026.
#' Habitat Affinity of Riverine Dissolved Organic Matter Linked to Molecular Traits.
#' *Global Change Biology*. [https://onlinelibrary.wiley.com/doi/10.1111/gcb.70736?af=R](https://onlinelibrary.wiley.com/doi/10.1111/gcb.70736?af=R); [PDF](https://jjwang.name/data/uploads/publications/Cui-et-al-2026-GCB.pdf)
#'
#' @examples 
#' \dontrun{
#' iHaf(data.ra.habi1, data.ra.habi2)
#' }
#' @rdname iHaf
#' @export 

iHaf <- function (data.ra.habi1, data.ra.habi2) {
  
  data_ra <- dplyr::bind_rows(data.ra.habi1, data.ra.habi2)
  data_ra[is.na(data_ra)] <- 0
  data_ra$group <- c(rep("habitat1", nrow(data.ra.habi1)), rep("habitat2", nrow(data.ra.habi2)))
  
  data_ra_habi1 <- subset(data_ra, group == "habitat1")
  data_ra_habi2 <- subset(data_ra, group == "habitat2")
  
  data_ra_habi1 <- data_ra_habi1 %>% dplyr::select(-c('group'))
  data_ra_habi2 <- data_ra_habi2 %>% dplyr::select(-c('group'))
  
  peak.list = colnames(data_ra)[1:(ncol(data_ra) - 1)]
  
  i = 1
  for (i in 1:length(peak.list)){
    
    dat.tmp = data.frame(
      peak = peak.list[i],
      m1 = mean(data_ra_habi1[, i], na.rm = T),
      sd1 = stats::sd(data_ra_habi1[, i], na.rm = T),
      n1 = length(data_ra_habi1[, i]),
      m2 = mean(data_ra_habi2[, i], na.rm = T),
      sd2 = stats::sd(data_ra_habi2[, i], na.rm = T),
      n2 = length(data_ra_habi2[, i])
    )
    
    iHaf.out.tmp = metafor::escalc(
      measure = "SMD",
      m1i = m1, sd1i = sd1, n1i = n1,
      m2i = m2, sd2i = sd2, n2i = n2,
      vtype = "UB",
      data = dat.tmp
    )
    iHaf.out.tmp <- iHaf.out.tmp[, c("peak", "yi")]
    
    iHaf.out.tmp$pvalue <- stats::kruskal.test(data_ra[,i] ~ group, data = data_ra)$p.value
    
    if (i == 1) {
      iHaf.out = iHaf.out.tmp
    } else {
      iHaf.out = rbind(iHaf.out, iHaf.out.tmp)
    }
  }
  colnames(iHaf.out)[2] = "iHaf"
  rownames(iHaf.out) <- iHaf.out$peak
  
  iHaf_wm.habi1 <- FD::functcomp(iHaf.out["iHaf"], as.matrix(data_ra_habi1))
  iHaf_wm.habi2 <- FD::functcomp(iHaf.out["iHaf"], as.matrix(data_ra_habi2))
  colnames(iHaf_wm.habi1) <- "iHaf_wm"
  colnames(iHaf_wm.habi2) <- "iHaf_wm"
  
  return(list(iHaf = iHaf.out, iHaf_wm_habi1 = iHaf_wm.habi1, iHaf_wm_habi2 = iHaf_wm.habi2))
  
}
