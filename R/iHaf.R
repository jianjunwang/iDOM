
#' @title The Indicator of Molecular Habitat Affinity (iHaf)
#' @description The iHaf function calculates the Indicator of Molecular Habitat Affinity (iHaf) based on molecular relative abundance data of habitat1 (data.ra.habi1) and habitat2 (data.ra.habi2).
#' @param data.ra.habi1 A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @param data.ra.habi2 A data frame containing the relative abundance of molecules, where rows represent samples and columns represent different molecules.
#' @return iHaf: A data frame of iHaf values for molecules.
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  iHaf(data.ra.habi1,data.ra.habi2)
#'  }
#' }
#' @rdname iHaf
#' @export 

iHaf <- function (data.ra.habi1, data.ra.habi2) 
{
  
  library(dplyr)
  library(metafor)
  
  data_ra <- bind_rows(data.ra.habi1, data.ra.habi2)
  data_ra[is.na(data_ra)] <- 0
  data_ra$group <- c(rep("habitat1", nrow(data.ra.habi1)), rep("habitat2", nrow(data.ra.habi2)))
  
  data_ra_habi1 <- subset(data_ra, group == "habitat1")
  data_ra_habi2 <- subset(data_ra, group == "habitat2")
  
  data_ra_habi1 <- data_ra_habi1 %>% select(-c('group'))
  data_ra_habi2 <- data_ra_habi2 %>% select(-c('group'))
  
  peak.list = colnames(data_ra)[1:(ncol(data_ra) - 1)]
  
  i = 1
  for (i in 1:length(peak.list)){

    dat.tmp = data.frame(
      peak = peak.list[i],
      m1 = mean(data_ra_habi1[, i], na.rm = T),
      sd1 = sd(data_ra_habi1[, i], na.rm = T),
      n1 = length(data_ra_habi1[, i]),
      m2 = mean(data_ra_habi2[, i], na.rm = T),
      sd2 = sd(data_ra_habi2[, i], na.rm = T),
      n2 = length(data_ra_habi2[, i])
    )
    
    iHaf.out.tmp = escalc(measure = "SMD",
                          m1i = m1, sd1i = sd1, n1i = n1,
                          m2i = m2, sd2i = sd2, n2i = n2,
                          vtype = "UB",
                          data = dat.tmp)
    iHaf.out.tmp <- iHaf.out.tmp[, c("peak", "yi")]
    
    iHaf.out.tmp$pvalue <- kruskal.test(data_ra[,i]~group, data = data_ra)$p.value
    
    if (i == 1) {
      iHaf.out = iHaf.out.tmp
    } else {
      iHaf.out = rbind(iHaf.out, iHaf.out.tmp)
    }
  }
  colnames(iHaf.out)[2] = "iHaf"
  
  return(iHaf.out)
  
}

