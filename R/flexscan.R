#' Analyze spatial count data using the flexible spatial scan statistic
#'
#' The rflexscan package provides functions and classes to analyze spatial count
#' data using the flexible spatial scan statistic developed by Tango and Takahashi (2005).
#' The original FleXScan software developed by Takahashi et al. is also 
#' available at \url{https://sites.google.com/site/flexscansoftware/home}.
#' 
#' @references 
#' \itemize{
#'   \item Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan statistic for detecting clusters, International Journal of Health Geographics 4:11.
#'   \item Takahashi K, Yokoyama T and Tango T. (2010). FleXScan v3.1: Software for the Flexible Scan Statistic. National Institute of Public Health, Japan.
#' }
#' 
"_PACKAGE"

#' @useDynLib rflexscan, .registration=TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' @export
flexscan.model <- c("POISSON", "BINOMIAL")

#' @export
flexscan.stattype <- c("ORIGINAL", "RESTRICTED")

#' @export
flexscan.scanmethod <- c("FLEXIBLE", "CIRCULAR")

#' @export
flexscan.rantype <- c("MULTINOMIAL", "POISSON")


#' Analyze spatial count data
#' 
#' Analyzes spatial count data using the flexible/circular 
#' spatial scan statistic.
#' 
#' @param case
#' The frequency of disease in each area.
#' The first column is the observed number of diseases, and the
#' second column is the expected number of diseases under the null
#' hypothesis, or the background population at risk in each area. 
#' 
#' @param coordinates
#' The geographic coordinates for each area.
#' Coordinates may be specified either using the standard Cartesian 
#' coordinates system or in latitude and longitude. The first column is the
#' x-coordinate or the latitude, and the second column is the y-coordinate or
#' longitude.
#' 
#' @param adj_mat
#' The adjacency matrix.
#' 
#' @param name
#' The name of each area.
#' 
#' @param clustersize
#' The number of maximum spatial cluster size to scan.
#' 
#' @param radius
#' Radius of Earth to calculate a distance between two sets of latitude and
#' longitude. It is approximately 6370 km in Japan.
#' 
#' @param model
#' Statistical model to be used (case-insensitive).
#' \describe{
#'   \item{"POISSON"}{for the data of the 'observed number' and the 'expected number'}
#'   \item{"BINOMIAL"}{for the data of the 'observed number' and the 'population'}
#' }
#' 
#' @param stattype
#' Statistic type to be used (case-insensitive).
#' \describe{
#'   \item{"ORIGINAL"}{the likelihood ratio statistic by Kulldorff and Nagarwalla (1995)}
#'   \item{"RESTRICTED"}{the restricted likelihood ratio statistic by Tango (2008), with a preset parameter \code{ralpha} for restriction}
#' }
#' 
#' @param scanmethod
#' Scanning method to be used (case-insensitive).
#' \describe{
#'   \item{"FLEXIBLE"}{flexible scan statistic by Tango and Takahashi (2005)}
#'   \item{"CIRCULAR"}{circular scan statistic by Kulldorff (1997)}
#' }
#' 
#' @param ralpha
#' Parameter for the restricted likelihood ratio statistic.
#' 
#' @param cartesian
#' Set \code{TRUE} for Cartesian coordinates.
#' 
#' @param simcount
#' The number of Monte Carlo replications to calculate a p-value for statistical
#' test.
#' 
#' @param rantype
#' The type of random number for Monte Carlo simulation (case-insensitive).
#' \describe{
#'   \item{"MULTINOMIAL"}{Total number of cases in whole area is fixed. It can be chosen in either Poisson or Binomial model.}
#'   \item{"POISSON"}{Total number of cases is not fixed. It can be chosen in Poisson model.}
#' }
#' 
#' @param ranseed
#' The seed for generating random numbers in the Monte Carlo simulation.
#' 
#' @param comments
#' Comments for the analysis which will be written in summary.
#' 
#' @return 
#' An \code{rflexscan} object which contains analysis results and specified
#' parameters.
#' 
#' @examples
#' # load sample data
#' data("saitama_heart_m")
#' data("saitama_coordinates")
#' data("saitama_adj_mat")
#' 
#' # Run FleXScan with default parameters
#' fls <- flexscan(saitama_heart_m, saitama_coordinates, saitama_adj_mat)
#' 
#' # Print summary to the terminal
#' summary(fls)
#' 
#' # Plot graph
#' plot(fls)
#' labs <- 1:length(fls$cluster)
#' legend("topright", legend = labs, col = rainbow(length(fls$cluster)), lty = 1)
#' 
#' @references
#' \itemize{
#'   \item Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan statistic for detecting clusters, International Journal of Health Geographics 4:11.
#'   \item Kulldorff M. and Nagarwalla N. (1995). Spatial disease clusters: Detection and Inference. Statistics in Medicine 14:799-810.
#'   \item Kulldorff M. (1997). A spatial scan statistic. Communications in Statistics: Theory and Methods, 26:1481-1496.
#'   \item Tango T. (2008). A spatial scan statistic with a restricted likelihood ratio. Japanese Journal of Biometrics 29(2):75-95.
#' }
#' 
#' @export
#' 
flexscan <- function(case,
                     coordinates,
                     adj_mat,
                     name=row.names(case),
                     clustersize=15,
                     radius=6370,
                     model=flexscan.model,
                     stattype=flexscan.stattype,
                     scanmethod=flexscan.scanmethod,
                     ralpha=0.2,
                     cartesian=FALSE,
                     simcount=999,
                     rantype=flexscan.rantype,
                     ranseed=4586111,
                     comments="",
                     verbose=FALSE) {
  model <- match.arg(toupper(model), flexscan.model)
  stattype <- match.arg(toupper(stattype), flexscan.stattype)
  scanmethod <- match.arg(toupper(scanmethod), flexscan.scanmethod)
  rantype <- match.arg(toupper(rantype), flexscan.rantype)
    
  row.names(case) <- name
  row.names(coordinates) <- name
  row.names(adj_mat) <- name
  
  if (cartesian) {
    colnames(coordinates) <- c("X", "Y")
  } else {
    colnames(coordinates) <- c("Latitude", "Longitude")
  }
  
  casefile <- tempfile()
  coofile <- tempfile()
  mt0file <- tempfile()
  mtrfile <- tempfile()
  resultfile <- tempfile()
  lambdafile <- tempfile()
  edgefile <- tempfile()
  nodefile <- tempfile()
  rfile <- tempfile()
  settingfile <- tempfile()
  stdoutfile <- tempfile()
  
  write.table(case, file = casefile, quote = FALSE, col.names = FALSE)
  write.table(coordinates, file = coofile, quote = FALSE, col.names = FALSE)
  
  diag(adj_mat) <- 2
  write.table(adj_mat, file = mt0file, quote = FALSE, col.names = FALSE)
  
  for (i in 1:nrow(adj_mat)) {
    cat(name[adj_mat[i,] != 0], "\n", file = mtrfile, append = TRUE)
  }
  
  # Write setting file
  cat("[FLEXSCAN]\n", sep = "", file = settingfile, append = TRUE)
  cat("VERSION=3.1.2\n", sep = "", file =settingfile, append = TRUE)
  cat("CASE=", casefile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("COORDINATES=", coofile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("MATRIX=", mt0file, "\n", sep = "", file =settingfile, append = TRUE)
  cat("MATRIXmtr=", mtrfile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RESULTS=", resultfile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("LAMBDAFILE=", lambdafile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("EDGEFILE=", edgefile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("NODEFILE=", nodefile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RFILE=", rfile, "\n", sep = "", file =settingfile, append = TRUE)
  cat("CLUSTERSIZE=", clustersize, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RADIUS=", radius, "\n", sep = "", file =settingfile, append = TRUE)
  cat("MODEL=", toupper(model), "\n", sep = "", file =settingfile, append = TRUE)
  cat("STATTYPE=", as.integer(stattype == "RESTRICTED"), "\n", sep = "", file =settingfile, append = TRUE)
  cat("SCANMETHOD=", toupper(scanmethod), "\n", sep = "", file =settingfile, append = TRUE)
  cat("RALPHA=", ralpha, "\n", sep = "", file =settingfile, append = TRUE)
  cat("CARTESIAN=", as.integer(cartesian), "\n", sep = "", file =settingfile, append = TRUE)
  cat("SIMCOUNT=", simcount, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RANTYPE=", toupper(rantype), "\n", sep = "", file =settingfile, append = TRUE)
  cat("RANSEED=", ranseed, "\n", sep = "", file =settingfile, append = TRUE)
  cat("COMMENT=", comments, "\n", sep = "", file =settingfile, append = TRUE)
  
  if (!verbose)
    sink(stdoutfile)
  
  start <- date()
  exit_code <- runFleXScan(settingfile)
  end <- date()
  
  if (!verbose)
    sink()
  
  result <- scan(resultfile, what = character(), sep = "\n", blank.lines.skip = FALSE)
  clst <- read.table(rfile, header = TRUE, stringsAsFactors = FALSE)
  
  if (toupper(model) == "POISSON") {
    clst <- apply(clst, 1, function(x){
      obj <- list()
      obj$max_dist <- as.numeric(x[1])
      obj$from <- as.character(x[2])
      obj$to <- as.character(x[3])
      obj$n_case <- as.integer(x[4])
      obj$expected <- as.numeric(x[5])
      obj$RR <- as.numeric(x[6])
      obj$stats <- as.numeric(x[7])
      obj$rank <- as.integer(x[8])
      obj$pval <- as.numeric(x[9])
      obj$area <- which(as.logical(as.integer(x[-(1:9)])))
      class(obj) <- 'Cluster'
      obj})
    colnames(case) <- c("Observed", "Expected")
  } else {
    clst <- apply(clst, 1, function(x){
      obj <- list()
      obj$max_dist <- as.numeric(x[1])
      obj$from <- as.character(x[2])
      obj$to <- as.character(x[3])
      obj$n_case <- as.integer(x[4])
      obj$population <- as.integer(x[5])
      obj$stats <- as.numeric(x[6])
      obj$rank <- as.integer(x[7])
      obj$pval <- as.numeric(x[8])
      obj$area <- which(as.logical(as.integer(x[-(1:8)])))
      class(obj) <- 'Cluster'
      obj})
    colnames(case) <- c("Observed", "Population")
  }
  
  diag(adj_mat) <- 0
  for (i in 1:length(clst)) {
    x <- clst[[i]]
    adj_mat[x$area,x$area] <- adj_mat[x$area,x$area] * (10 * i)
  }

  retval <- list(case=case, coordinates=coordinates, name=name,
                 cluster=clst, clustersize=clustersize, 
                 radius=radius, model=model, stattype=stattype,
                 scanmethod=scanmethod, ralpha=ralpha, cartesian=cartesian,
                 simcount=simcount, rantype=rantype, ranseed=ranseed,
                 comments=comments, summary=result, adj_mat=adj_mat)
  class(retval) <- "rflexscan"

  return(retval)
}

#' Print summary to the terminal
#' 
#' @export
#' 
summary.rflexscan <- function(object, ...) {
  cat(object$summary, sep = "\n")
}

#' Graph plotting of flexscan results
#' 
#' @examples
#' # highlight all clusters
#' plot(flexscan.output)
#' 
#' # highlight clusters with rank 1 and 2
#' plot(flexscan.output, rank = c(1,2))
#' 
#' # highlight clusters of P-value <= 0.05
#' plot(flexscan.output, pval = 0.05)
#' 
#' @import igraph
#' 
#' @export
#' 
plot.rflexscan <- function(object,
                           xlim=c(min(object$coordinates[,1]), max(object$coordinates[,1])),
                           ylim=c(min(object$coordinates[,2]), max(object$coordinates[,2])),
                           col=rainbow(length(object$cluster)),
                           rank=1:length(object$cluster),
                           pval=1,
                           xlab=colnames(object$coordinates)[1],
                           ylab=colnames(object$coordinates)[2],
                           ...) {
  g <- graph_from_adjacency_matrix(object$adj_mat, mode = "undirected", diag = FALSE, weighted = TRUE)
  V(g)$size <- 1
  V(g)$frame.color <- "gray40"
  V(g)$color <- "white"
  V(g)$label <- ""
  E(g)$color <- "gray40"
  
  labs <- character()
  # color clusters
  for (i in 1:length(object$cluster)) {
    if (i %in% rank & object$cluster[[i]]$pval <= pval) {
      V(g)$color[object$cluster[[i]]$area] <- col[i]
      E(g)$color[E(g)$weight == 10 * i] <- col[i]
      labs <- c(labs, sprintf("%d (P=%f)", i, object$cluster[[i]]$pval))
    }
  }
  
  if (object$cartesian) {
    plot(g, axes = TRUE, layout = as.matrix(object$coordinates[,c(1,2)]), rescale = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  } else {
    # flip X-Y (x:longitude, y:latitude)
    plot(g, axes = TRUE, layout = as.matrix(object$coordinates[,c(2,1)]), rescale = FALSE, xlab = ylab, ylab = xlab, xlim = ylim, ylim = xlim, ...)
  }
}
