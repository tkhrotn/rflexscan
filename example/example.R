library(rflexscan)
library(spdep)
data("nc.sids")
 
# calculate the expected numbers of cases
expected <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)
 
# run FleXScan with default parameters
fls <- flexscan(lat = nc.sids$lat, lon = nc.sids$lon,
                 observed = nc.sids$SID74,
                 expected = expected,
                 name = rownames(nc.sids),
                 nb = ncCR85.nb, verbose = TRUE)
