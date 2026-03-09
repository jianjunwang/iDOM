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

load("data-raw/mol.data.ra_habitat.Rdata")
dim(data.ra.habi1)
mol.data.ra.habi1 <- data.ra.habi1

dim(data.ra.habi2)
mol.data.ra.habi2 <- data.ra.habi2

usethis::use_data(mol.data.ra.habi1, overwrite = TRUE)
usethis::use_data(mol.data.ra.habi2, overwrite = TRUE)

# write data in correct format to data folder ----
usethis::use_data(mol.data, overwrite = TRUE)
usethis::use_data(mol.trait, overwrite = TRUE)
usethis::use_data(envi, overwrite = TRUE)
usethis::use_data(mol.dark.matter, overwrite = TRUE)
usethis::use_data(micro.data, overwrite = TRUE)
usethis::use_data(Transformation_Database, overwrite = TRUE)
