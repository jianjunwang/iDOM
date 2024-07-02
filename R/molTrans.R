
#' @title Molecular Transformation Analysis
#' @description This function analyzes molecular data to identify potential transformations between molecular peaks based on their mass differences. 
#' It supports both dataset-level and sample-level analysis, enabling the examination of transformations within a comprehensive dataset or within individual samples respectively.
#' @param mol.data Data frame containing molecular data where rows represent different samples or conditions, and columns represent individual molecules. Default: mol.data.
#' @param mol.trait Data frame containing molecular trait data such as mass or other properties associated with each molecule. Default: mol.trait.
#' @param error.term Numeric value representing the allowable error margin in the mass difference calculations. Default: 1e-05.
#' @param type Character string specifying the level of analysis; accepts "Dataset" for a comprehensive analysis across all samples, or "Sample" for an analysis focused on individual samples. Default: 'Dataset'.
#' @return A list containing two elements: `Peak.2.peak`, a data frame of peak-to-peak relationships highlighting the transformations, and `Peak.profile.dataset` or `Peak.profile.sample` depending on the analysis type, which summarizes the number of transformations each peak is involved in.
#' @rdname molTrans
#' @export

molTrans <- function(mol.data, mol.trait, error.term = 0.000010, type = "Dataset"){
  
  library(tidyverse)
  
  options(digits = 10) # Significant figures in mass resolution data
  
  # Checking row names consistency
  if(identical(x = colnames(mol.data), y = rownames(mol.trait)) == FALSE){
    stop("Something is incorrect in your row names")
  }
  
  # Loading in transformations
  trans.full = Transformation_Database
  
  if (type == "Dataset") {
    
    # Creating a presence matrix for all peaks observed within the dataset
    bulk.peaks = as.data.frame(cbind("Sample_All",row.names(mol.trait)))
    bulk.peaks[,2] = as.numeric(bulk.peaks[,2])
    colnames(bulk.peaks) = c("sample", "peak.x")
    
    # Sort molecules from smallest to largest by mass
    bulk.peaks <- bulk.peaks[order(bulk.peaks$peak.x),]
    
    peak.2.peak = NULL
    
    i = 1
    # Running a loop to compare each peak to each other peak
    for(i in 1:(ncol(mol.data)-1)){ # I cannot stress the importance of the "-1" here...
      
      # Creating a data matrix to ensur no repeat or negative differences
      Distance_Results = bulk.peaks[-1:-i,] # Removing all peaks up to, and including the current peak
      Distance_Results$peak.y = bulk.peaks$peak.x[i] # Setting the peak of interest
      Distance_Results$Dist = Distance_Results$peak.x - Distance_Results$peak.y # Finding the difference between all peaks and the peak of interest
      
      # Adding in error terms to the matrix
      Distance_Results$Dist.plus = Distance_Results$Dist + error.term
      Distance_Results$Dist.minus = Distance_Results$Dist - error.term
      Distance_Results$Trans.name = -999
      
      for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
        mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
        if (length(mass.diff) > 1) { break() }
        Distance_Results$Trans.name[which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)] = current.trans
      }
      
      # Removing differences that didn't match any transformation
      Distance_Results = Distance_Results %>% 
        filter(!Trans.name == -999)
      head(Distance_Results)
      
      # Building a larger peak.2.peak file
      peak.2.peak = rbind(peak.2.peak, Distance_Results)
      
      print(paste("Finished running through peak #", i, " on ", date(), sep = ""))
      
    }
    
    # Creating a num.trans file for network generation
    peak.stack = as.data.frame(c(peak.2.peak$peak.x, peak.2.peak$peak.y)); head(peak.stack)
    peak.profile = as.data.frame(tapply(X = peak.stack[,1], INDEX = peak.stack[,1], FUN = 'length' )); dim(peak.profile)
    colnames(peak.profile) = 'num.trans.involved.in'
    peak.profile$peak = row.names(peak.profile)
    peak.profile$sample = "Sample_All"
    
    # "num.trans.involved.in" is zero 
    mol.trans0 = mol.trait[!(rownames(mol.trait) %in% rownames(peak.profile)),]
    peak.profile.trans0 = data.frame(num.trans.involved.in = 0, peak = rownames(mol.trans0),sample = "Sample_All")
    rownames(peak.profile.trans0) = peak.profile.trans0$peak
    
    peak.profile.dataset = rbind(peak.profile, peak.profile.trans0)
    peak.profile.dataset = peak.profile.dataset[rownames(mol.trait),]
    
    Transformation.results <- list(Peak.2.peak = peak.2.peak, Peak.profile.dataset = peak.profile.dataset)
    return(Transformation.results)
    
  }else if (type == "Sample"){
    
    peak.2.peak = NULL
    peak.profile.sample = NULL
    
    # pull out just the sample names
    samples.to.process = rownames(mol.data)
    
    # current.sample = samples.to.process[1]
    
    for (current.sample in samples.to.process) {
      
      one.sample.matrix = as.data.frame(t(mol.data[which(rownames(mol.data) == current.sample),,drop = F]))
      one.sample.matrix = data.frame(peak = rownames(one.sample.matrix),one.sample.matrix)
      
      Sample_Peak_Mat <- one.sample.matrix %>% gather("sample", "value", -1) %>% filter(value > 0) %>% select(sample, peak)
      
      Distance_Results <- Sample_Peak_Mat %>% 
        left_join(Sample_Peak_Mat, by = "sample",relationship = "many-to-many") %>% 
        mutate(peak.x = as.numeric(peak.x),peak.y = as.numeric(peak.y)) %>% 
        filter(peak.x > peak.y) %>% mutate(Dist = peak.x - peak.y) %>% 
        select(sample,peak.x,peak.y,Dist)
      
      Distance_Results$Dist.plus = Distance_Results$Dist + error.term
      Distance_Results$Dist.minus = Distance_Results$Dist - error.term
      Distance_Results$Trans.name = -999
      head(Distance_Results)
      
      dist.unique = unique(Distance_Results[,'sample']) #unique samples
      
      for (current.trans in unique(trans.full$Name)) { # note that for masses with multiple names, only the last name is going to be recorded
        mass.diff = trans.full$Mass[which(trans.full$Name == current.trans)]
        if (length(mass.diff) > 1) { break() }
        Distance_Results$Trans.name[which(Distance_Results$Dist.plus >= mass.diff & Distance_Results$Dist.minus <= mass.diff)] = current.trans
      }
      
      Distance_Results = Distance_Results %>% 
        filter(!Trans.name == -999)
      head(Distance_Results)
      
      # find the number of transformations each peak was associated with
      peak.stack = as.data.frame(c(Distance_Results$peak.x,Distance_Results$peak.y)); head(peak.stack)
      peak.profile = as.data.frame(tapply(X = peak.stack[,1],INDEX = peak.stack[,1],FUN = 'length' )); dim(peak.profile)
      colnames(peak.profile) = 'num.trans.involved.in'
      peak.profile$peak = row.names(peak.profile)
      peak.profile$sample = dist.unique
      head(peak.profile)
      
      # "num.trans.involved.in" is zero 
      mol.trans0 = Sample_Peak_Mat[!(Sample_Peak_Mat$peak %in% rownames(peak.profile)),]
      peak.profile.trans0 = data.frame(num.trans.involved.in = 0, peak = mol.trans0$peak, sample = dist.unique)
      rownames(peak.profile.trans0) = peak.profile.trans0$peak
      
      peak.profile.all = rbind(peak.profile, peak.profile.trans0)
      peak.profile.all = peak.profile.all[Sample_Peak_Mat$peak,];dim(peak.profile.all)
      
      peak.2.peak = rbind(peak.2.peak,Distance_Results)
      peak.profile.sample = rbind(peak.profile.sample, peak.profile.all)
      
      print(paste("Finished running through sample #", dist.unique, " on ", date(), sep = ""))
    }
    
    Transformation.results <- list(Peak.2.peak = peak.2.peak, Peak.profile.sample = peak.profile.sample)
    return(Transformation.results)

  }
}

