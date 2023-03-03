## code to prepare `DATASET` dataset goes here

usethis::use_data(DATASET, overwrite = TRUE)

rm(list = ls())
ls()

my_file <- read.csv("R/my_data.csv",header = T)
ls()

csv = system.file("extdata","cdr3.csv", package = "ClusTCR2")
read.csv(csv)
