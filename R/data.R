
#' Example datasets included in iDOM
#'
#' The `iDOM` package provides example datasets of dissolved organic matter (DOM) and microbial communities under experimental warming.
#' These datasets were generated from a laboratory microcosm experiment using sterilized Taihu Lake sediments as the organic carbon source, inoculated with distinct microbial communities from sediments of subtropical and temperate lakes in China.
#' The microcosms were incubated in the dark for 1 month at six temperature levels (5, 10, 15, 20, 25, and 30 degrees C), with three replicates per temperature treatment, resulting in 36 samples across two climate zones.
#' 
#' Five core example datasets are included: molecular intensity data, molecular element composition data, environmental data, dark matter data, and microbial data. 
#' In addition, the package includes a molecular transformation database and two habitat-specific relative abundance datasets for examples involving habitat affinity analysis.
#'
#' @format `mol.data`: A data frame with 36 samples in rows and 5474 assigned DOM molecules in columns, giving molecular intensity data.
#' @format `mol.trait`: A data frame with 5474 molecules in rows and 29 columns containing molecular element composition assignments and derived molecular traits for the molecules in `mol.data`.
#' @format `envi`: A data frame with 36 samples in rows and 3 columns describing sample metadata for `mol.data`. Variables are `Source`, `Temperature`, and `Parallel`.
#' @format `mol.dark.matter`: A data frame with 36 samples in rows and 5779 unassigned molecular formula in columns, representing molecular dark matter intensity data.
#' @format `micro.data`: A data frame with 36 samples in rows and 463 microbial bacterial genera in columns, giving microbial relative abundance data.
#' @format `Transformation_Database`: A data frame with 1255 rows and 2 columns, `Name` and `Mass`, providing transformation names and their corresponding mass differences.
#' @format `mol.data.ra.habi1`: A data frame with 50 samples in rows and 1000 molecules in columns representing a habitat-specific DOM relative abundance matrix.
#' @format `mol.data.ra.habi2`: A data frame with 50 samples in rows and 1000 molecules in columns representing a second habitat-specific DOM relative abundance matrix.
#'
#' @details
#' The row names of `mol.data`, `envi`, `mol.dark.matter`, and `micro.data` identify the same set of 36 experimental samples.
#' The row names of `mol.trait` correspond to the molecular mass labels used as column names in `mol.data`. 
#' Together, the five core datasets can be used to demonstrate analyses of DOM composition, molecular traits, environmental responses, dark matter effects, and DOM-microbe associations. 
#' 
#' The objects `mol.data.ra.habi1` and `mol.data.ra.habi2` are included for examples involving habitat affinity analysis.
#'
#' @source Internal example datasets distributed with the `iDOM` package.
#' @name iDOM-data
#' @aliases mol.data mol.trait envi mol.dark.matter micro.data Transformation_Database mol.data.ra.habi1 mol.data.ra.habi2
NULL
