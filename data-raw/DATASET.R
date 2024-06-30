## code to prepare `DATASET` dataset goes here

# usethis::use_data(DATASET, overwrite = TRUE)

rm(list = ls())
load("data-raw/mol.data.RData")
dim(mol.data)

load("data-raw/mol.trait.RData")
dim(mol.trait)

load("data-raw/envi.RData")
dim(envi)

load("data-raw/mol.dark.matter.RData")
dim(mol.dark.matter)

load("data-raw/micro.data.RData")
dim(micro.data)

load("data-raw/Transformation_Database.RData")
dim(Transformation_Database)

# write data in correct format to data folder ----
usethis::use_data(mol.data, overwrite = TRUE)
usethis::use_data(mol.trait, overwrite = TRUE)
usethis::use_data(envi, overwrite = TRUE)
usethis::use_data(mol.dark.matter, overwrite = TRUE)
usethis::use_data(micro.data, overwrite = TRUE)
usethis::use_data(Transformation_Database, overwrite = TRUE)
