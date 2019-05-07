#' Analyze spatial count data using the flexible spatial scan statistic
#'
#' The rflexscan package provides functions and classes to analyze spatial count
#' data using the flexible spatial scan statistic developed by Tango and 
#' Takahashi (2005). This package designed for any of the following interrelated
#'  purposes:
#' \enumerate{
#'   \item To evaluate reported spatial disease clusters, to see if they are 
#'         statistically significant.
#'   \item To test whether a disease is randomly distributed over space.
#'   \item To perform geographical surveillance of disease, to detect areas of 
#'         significantly high rates.
#' }
#' 
#' @references 
#' \itemize{
#'   \item Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan 
#'   statistic for detecting clusters, International Journal of Health
#'   Geographics 4:11.
#'   \item Takahashi K, Yokoyama T and Tango T. (2010). FleXScan v3.1: Software 
#'   for the Flexible Scan Statistic. National Institute of Public Health, Japan,
#'   \url{https://sites.google.com/site/flexscansoftware/home}.
#' }
#' 
#' @seealso \code{\link{flexscan}}
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


#' Detect spatial disease clusters using the flexible/circular scan statistic
#' 
#' This function Analyzes spatial count data using the flexible spatial scan 
#' statistic developed by Tango and Takahashi (2005) and Kulldorff‟s circular 
#' spatial scan statistic (1997) and detect spatial disease clusters.
#' 
#' @param x
#' An array of X-coordinates or a column name in \code{data}.
#' 
#' @param y
#' An arrya of Y-coordinates or a column name in \code{data}.
#' 
#' @param lat
#' An array of latitude or a column name in \code{data}.
#' 
#' @param lon
#' An array of longitude or a column name in \code{data}.
#' 
#' @param observed
#' An array of observed number of diseases or a clumn name in \code{data}.
#' 
#' @param expected
#' An array of expected number of diseases under the null hypothesis or a clumn
#' name in `data`. This is used on "Poisson" model.
#' 
#' @param population
#' An array of background population at risk in each area or a clumn name in 
#' `data`. This is used on "Binomial" model.
#' 
#' @param data
#' A dataset containing coordinates, names for each area, and disease case data.
#' 
#' @param nb
#' A neighbours list or an adjacency matrix.
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
#' @param stattype
#' Statistic type to be used (case-insensitive).
#' \describe{
#'   \item{"ORIGINAL"}{the likelihood ratio statistic by Kulldorff and
#'   Nagarwalla (1995)}
#'   \item{"RESTRICTED"}{the restricted likelihood ratio statistic by Tango 
#'   (2008), with a preset parameter \code{ralpha} for restriction}
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
#' @param simcount
#' The number of Monte Carlo replications to calculate a p-value for statistical
#' test.
#' 
#' @param rantype
#' The type of random number for Monte Carlo simulation (case-insensitive).
#' \describe{
#'   \item{"MULTINOMIAL"}{Total number of cases in whole area is fixed. It can 
#'   be chosen in either Poisson or Binomial model.}
#'   \item{"POISSON"}{Total number of cases is not fixed. It can be chosen in 
#'   Poisson model.}
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
#'   Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan 
#'   statistic for detecting clusters, International Journal of Health 
#'   Geographics 4:11.
#'   
#'   Kulldorff M. and Nagarwalla N. (1995). Spatial disease clusters: 
#'   Detection and Inference. Statistics in Medicine 14:799-810.
#'   
#'   Kulldorff M. (1997). A spatial scan statistic. Communications in 
#'   Statistics: Theory and Methods, 26:1481-1496.
#'   
#'   Tango T. (2008). A spatial scan statistic with a restricted 
#'   likelihood ratio. Japanese Journal of Biometrics 29(2):75-95.
#' }
#' 
#' @seealso \code{\link{summary.rflexscan}}, \code{\link{plot.rflexscan}}
#' 
#' @export
#' 
flexscan <- function(x, y, lat, lon, nb,
                     name, observed, expected, population, data,
                     clustersize=15,
                     radius=6370,
                     stattype=flexscan.stattype,
                     scanmethod=flexscan.scanmethod,
                     ralpha=0.2,
                     simcount=999,
                     rantype=flexscan.rantype,
                     ranseed=4586111,
                     comments="",
                     verbose=FALSE) {
  call <- match.call()

  stattype <- match.arg(toupper(stattype), flexscan.stattype)
  scanmethod <- match.arg(toupper(scanmethod), flexscan.scanmethod)
  rantype <- match.arg(toupper(rantype), flexscan.rantype)

  if (missing(data)) {
    if (!missing(lat) && !missing(lon)) {
      coordinates <- cbind(lat, lon)
      latlon <- TRUE
    } else {
      coordinates <- cbind(x, y)
      latlon <- FALSE
    }
    
    if (!missing(expected)) {
      case <- cbind(observed, expected)
      model <- "POISSON"
    } else {
      case <- cbind(observed, population)
      model <- "BINOMIAL"
    }
    
    row.names(coordinates) <- as.character(name)
    row.names(case) <- as.character(name)
  } else {
    if (!is.data.frame(data)) {
      data <- data.frame(data)
    }
    
    if (!missing(lat) && !missing(lon)) {
      coordinates <- cbind(data[,lat], data[,lon])
      latlon <- TRUE
    } else {
      coordinates <- cbind(data[,x], data[,y])
      latlon <- FALSE
    }
    
    if (!missing(expected)) {
      case <- cbind(data[,observed], data[,expected])
      model <- "POISSON"
    } else {
      case <- cbind(data[,observed], data[,population])
      model <- "BINOMIAL"
    }
    
    row.names(coordinates) <- as.character(data[,name])
    row.names(case) <- as.character(data[,name])
    name <- as.character(data[,name])
  } 
  
  if (is.matrix(nb)) {
    adj_mat <- nb
  } else {
    adj_mat <- matrix(0, nrow = nrow(coordinates), ncol = nrow(coordinates))
    for (i in 1:nrow(coordinates)) {
      adj_mat[i, nb[[i]]] <- 1
    }
  }
  row.names(adj_mat) <- row.names(coordinates)
  colnames(adj_mat) <- row.names(coordinates)
  
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
  cat("MODEL=", model, "\n", sep = "", file =settingfile, append = TRUE)
  cat("STATTYPE=", as.integer(stattype == "RESTRICTED"), "\n", sep = "", file =settingfile, append = TRUE)
  cat("SCANMETHOD=", scanmethod, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RALPHA=", ralpha, "\n", sep = "", file =settingfile, append = TRUE)
  cat("CARTESIAN=", as.integer(!latlon), "\n", sep = "", file =settingfile, append = TRUE)
  cat("SIMCOUNT=", simcount, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RANTYPE=", rantype, "\n", sep = "", file =settingfile, append = TRUE)
  cat("RANSEED=", ranseed, "\n", sep = "", file =settingfile, append = TRUE)
  cat("COMMENT=", comments, "\n", sep = "", file =settingfile, append = TRUE)
  
  if (!verbose)
    sink(stdoutfile)
  
  start <- date()
  exit_code <- runFleXScan(settingfile)
  end <- date()
  
  if (!verbose)
    sink()
  
  result <- scan(resultfile, what = character(), sep = "\n", blank.lines.skip = FALSE, quiet = TRUE)
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
  
  clusterrank <- numeric(nrow(case))
  for (i in 1:length(clst)) {
    clusterrank[clst[[i]]$area] <- i
  }
  
  diag(adj_mat) <- 0
  for (i in 1:length(clst)) {
    x <- clst[[i]]
    adj_mat[x$area,x$area] <- adj_mat[x$area,x$area] * (10 * i)
  }

  retval <- list(call=call, case=case, coordinates=coordinates, name=name,
                 cluster=clst, clusterrank=clusterrank, clustersize=clustersize, 
                 radius=radius, model=model, stattype=stattype,
                 scanmethod=scanmethod, ralpha=ralpha, latlon=latlon,
                 simcount=simcount, rantype=rantype, ranseed=ranseed,
                 comments=comments, summary=result, adj_mat=adj_mat)
  class(retval) <- "rflexscan"

  return(retval)
}


#' Summarizing rflexscan results
#' 
#' Print and summary method functions for rflexscan objects.
#' 
#' @export
#' 
summary.rflexscan <- function(object, ...) {
  n_cluster <- length(object$cluster)
  total_areas <- nrow(object$case)
  total_cases <- sum(object$case[,"Observed"])
  areas <- lapply(object$cluster, function(i) {i$area})

  n_area <- sapply(object$cluster, function(i){length(i$area)})
  max_dist <- sapply(object$cluster, function(i) {i$max_dist})
  n_case <- sapply(object$cluster, function(i) {i$n_case})
  stats <- sapply(object$cluster, function(i) {i$stats})
  pval <- sapply(object$cluster, function(i) {i$pval})
  
  if (toupper(object$model) == "POISSON") {
    expected <- sapply(object$cluster, function(i) {i$expected})
    RR <- sapply(object$cluster, function(i) {i$RR})

    table <- cbind(NumArea=n_area, MaxDist=max_dist, Case=n_case, 
                   Expected=expected, RR=RR, Stats=stats, P=pval)
  } else {
    population <- sapply(object$cluster, function(i) {i$population})

    table <- cbind(NumArea=n_area, MaxDist=max_dist, Case=n_case,
                   Population=population, Stats=stats, P=pval)
  }
  row.names(table) <- 1:n_cluster

  retval <- list(call=object$call, clustersize=object$clustersize, 
                 n_cluster=n_cluster, name=object$name,
                 total_areas=total_areas, total_cases=total_cases,
                 areas=areas, stattype=object$stattype, model=object$mode,
                 scanmethod=object$scanmethod, latlon=object$latlon,
                 clusters=table)
  
  class(retval) <- "summary.rflexscan"
  return(retval)
}



#' Print results of flexscan
#' 
#' @export
#' 
print.summary.rflexscan <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Clusters:\n")
  signif <- symnum(x$clusters[,"P"], corr = FALSE, 
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " "))
  table <- cbind(x$clusters, signif)
  colnames(table)[ncol(table)] <- ""
  print(table, quote = FALSE, right = TRUE, print.gap = 2)
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
  
  cat("Census areas in cluster:\n")
  for (i in 1:x$n_cluster) {
    cat(sprintf("[%d]", i), x$name[x$areas[[i]]], sep = " ", fill = TRUE)
  }
  
  cat("\nLimit length of cluster:", x$clustersize, "\n")
  cat("Number of census areas:", x$total_areas, "\n")
  cat("Total cases:", x$total_cases, "\n")
  if (x$latlon) {
    cat("Coordinates: Latitude/Longitude\n")
  } else {
    cat("Coordinates: Cartesian\n")
  }
  cat("Scanning Method:", x$scanmethod, "\n")
  cat("Model:", x$model, "\n\n")
}
  

#' Graph plotting of flexscan results
#' 
#' The \code{plot} method for flexscan objects.
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
                           col=rainbow(length(object$cluster)),
                           rank=1:length(object$cluster),
                           pval=1,
                           vertexsize=max(object$coordinates[,1])-min(object$coordinates[,1]),
                           xlab=colnames(object$coordinates)[1],
                           ylab=colnames(object$coordinates)[2],
                           xlim=c(min(object$coordinates[,1]), max(object$coordinates[,1])),
                           ylim=c(min(object$coordinates[,2]), max(object$coordinates[,2])),
                           ...) {
  g <- graph_from_adjacency_matrix(object$adj_mat, mode = "undirected", diag = FALSE, weighted = TRUE)
  V(g)$size <- vertexsize
  V(g)$frame.color <- "gray40"
  V(g)$color <- "white"
  V(g)$label <- ""
  E(g)$color <- "gray40"
  
  # color clusters
  for (i in 1:length(object$cluster)) {
    if (i %in% rank & object$cluster[[i]]$pval <= pval) {
      V(g)$color[object$cluster[[i]]$area] <- col[i]
      E(g)$color[E(g)$weight == 10 * i] <- col[i]
    }
  }
  
  if (!object$latlon) {
    plot(g, axes = TRUE, layout = as.matrix(object$coordinates[,c(1,2)]), rescale = FALSE, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  } else {
    # flip X-Y (x:longitude, y:latitude)
    plot(g, axes = TRUE, layout = as.matrix(object$coordinates[,c(2,1)]), rescale = FALSE, xlab = ylab, ylab = xlab, xlim = ylim, ylim = xlim, ...)
  }
}
