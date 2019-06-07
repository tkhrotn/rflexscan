saitama_coordinates <- read.table("C:/Users/otani/Desktop/Work/データ/H11_47ken-coo/11埼玉県.coo", row.names=1)
colnames(saitama_coordinates) <- c("Latitude", "Longitude")

mtr <- strsplit(scan("C:/Users/otani/Desktop/Work/データ/H11_47ken-mtr/11埼玉県.mtr", what = character(), sep = "\n", blank.lines.skip = FALSE), " ")

mtr2adjmat <- function(mtr, names) {
  adjmat <- sapply(mtr, function(x) {
    as.integer(names %in% x)
  })
  diag(adjmat) <- 0
  row.names(adjmat) <- names
  colnames(adjmat) <- names
  
  adjmat
}

saitama_adj_mat <- mtr2adjmat(mtr, row.names(saitama_coordinates))

data("saitama_heart_m")
row.names(saitama_heart_m) <- row.names(saitama_coordinates)

devtools::use_data(saitama_coordinates, saitama_adj_mat, saitama_heart_m, overwrite = TRUE)

data("saitama_coordinates")
names <- row.names(saitama_coordinates)
library(rgdal)
saitama_coordinates <- project(as.matrix(saitama_coordinates[,c(2,1)]), "+proj=utm +zone=54 ellps=WGS84")
colnames(saitama_coordinates) <- c("X", "Y")
row.names(saitama_coordinates) <- names