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
#' This package implements a wrapper for the C routine used in the FleXScan 3.1.2 
#' developed by Takahashi, Yokoyama, and Tango.
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
#' @seealso \code{\link{rflexscan}}
#' @aliases NULL rflexscan-package
#' 
"_PACKAGE"

#' @useDynLib rflexscan
#' @importFrom Rcpp sourceCpp
NULL

flexscan.model <- c("POISSON", "BINOMIAL")
flexscan.stattype <- c("ORIGINAL", "RESTRICTED")
flexscan.scanmethod <- c("FLEXIBLE", "CIRCULAR")
flexscan.rantype <- c("MULTINOMIAL", "POISSON")


#' Detect spatial disease clusters using the flexible/circular scan statistic
#' 
#' This function analyzes spatial count data using the flexible spatial scan 
#' statistic developed by Tango and Takahashi (2005) or Kulldorff's circular 
#' spatial scan statistic (1997), and detect spatial disease clusters.
#' 
#' @details
#' Centroid coordinates for each region should be specified EITHER by Cartesian 
#' coordinates using arguments \code{x} and \code{y} or by latitudes and 
#' longitudes using arguments \code{lat} and \code{lon}.
#' 
#' @param x
#' A vector of X-coordinates.
#' 
#' @param y
#' A vector of Y-coordinates.
#' 
#' @param lat
#' A vector of latitude.
#' 
#' @param lon
#' A vector of longitude.
#' 
#' @param observed
#' A vector with the observed number of disease cases.
#' 
#' @param expected
#' A vector with the expected number of disease cases under the null hypothesis. 
#' This is used on "Poisson" model.
#' 
#' @param population
#' A vector with the background population at risk in each area. 
#' This is used on "Binomial" model.
#' 
#' @param nb
#' A neighbors list or an adjacency matrix.
#' 
#' @param name
#' A vector of names of each area.
#' 
#' @param clustersize
#' The number of maximum spatial cluster size to scan, i.e., the maximum number 
#' of regions included in the detected cluster
#' 
#' @param radius
#' Radius of Earth to calculate a distance between two sets of latitude and
#' longitude. It is approximately 6370 km in Japan. This parameter is used when 
#' \code{lat} and \code{lon} are specified. This is DEPRECATED. The 
#' distance calculated using this parameter is not accurate. This feature is 
#' implemented to maintain compatibility with FleXScan. It is recommended to 
#' transform latitude and longitude onto the Cartesian coordinate system 
#' beforehand and use the x and y parameters.
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
#' Threshold parameter of the middle p-value for the restricted likelihood ratio
#' statistic.
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
#' @param comments
#' Comments for the analysis which will be written in summary.
#' 
#' @param verbose
#' Print progress messages.
#' 
#' @return 
#' An \code{rflexscan} object which contains analysis results and specified
#' parameters.
#' 
#' @examples
#' # load sample data (North Carolina SIDS data)
#' library(spdep)
#' data("nc.sids")
#' 
#' # calculate the expected numbers of cases
#' expected <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)
#' 
#' # run FleXScan
#' fls <- rflexscan(x = nc.sids$x, y = nc.sids$y,
#'                  observed = nc.sids$SID74,
#'                  expected = expected,
#'                  name = rownames(nc.sids),
#'                  clustersize = 10,
#'                  nb = ncCR85.nb)
#' 
#' # print rflexscan object
#' print(fls)
#' 
#' # print properties of the most likely cluster
#' print(fls$cluster[[1]])
#' 
#' # print summary to the terminal
#' summary(fls)
#' 
#' # plot graph
#' plot(fls, col = palette())
#' labs <- 1:length(fls$cluster)
#' legend("bottomleft", legend = labs, col = palette(), lty = 1)
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
#' 
#' @seealso \link{summary.rflexscan}, \link{plot.rflexscan}, \link{choropleth}
#' 
#' @export
#' 
rflexscan <- function(x, y, lat, lon, 
                      name, observed, expected, population, nb,
                      clustersize=15,
                      radius=6370,
                      stattype="ORIGINAL",
                      scanmethod="FLEXIBLE",
                      ralpha=0.2,
                      simcount=999,
                      rantype="MULTINOMIAL",
                      comments="",
                      verbose=FALSE) {
  call <- match.call()

  stattype <- match.arg(toupper(stattype), flexscan.stattype)
  scanmethod <- match.arg(toupper(scanmethod), flexscan.scanmethod)
  rantype <- match.arg(toupper(rantype), flexscan.rantype)
  
  # replace space
  name <- sub(" ", "_", name)

  if (!missing(lat) && !missing(lon) && missing(x) && missing(y)) {
    coordinates <- cbind(lat, lon)
    latlon <- TRUE
  } else if (missing(lat) && missing(lon) && !missing(x) && !missing(y)) {
    coordinates <- cbind(x, y)
    latlon <- FALSE
  } else {
    stop("Coordinates are not properly specified.")
  }
  
  if (missing(observed)) {
    stop("Observed numbers of diseases are not specified.")
  }
  
  if (!missing(expected)) {
    case <- cbind(observed, expected)
    model <- "POISSON"
  } else if (!missing(population)) {
    case <- cbind(observed, population)
    model <- "BINOMIAL"
  } else {
    stop("Expected numbers of diseases or background population are not specified.")
  }
  
  row.names(coordinates) <- as.character(name)
  row.names(case) <- as.character(name)

  if (missing(nb)) {
    stop("A neighbours list or an adjacency matrix are not specified.")
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
  diag(adj_mat) <- 2

  setting <- list()
  setting$clustersize <- clustersize
  setting$radius <- radius
  setting$model <- as.integer(model == "BINOMIAL")
  setting$stattype <- as.integer(stattype == "RESTRICTED")
  setting$scanmethod <- as.integer(scanmethod == "CIRCULAR")
  setting$ralpha <- ralpha
  setting$cartesian <- as.integer(!latlon)
  setting$simcount <- simcount
  setting$rantype <- as.integer(rantype == "POISSON")
      
  if (!verbose) {
    output <- capture.output({
      start <- date()
      clst <- runFleXScan(setting, case, coordinates, adj_mat)
      end <- date()
    })
  } else {
    start <- date()
    clst <- runFleXScan(setting, case, coordinates, adj_mat)
    end <- date()
  }
  
  if (toupper(model) == "POISSON") {
    colnames(case) <- c("Observed", "Expected")
  } else {
    colnames(case) <- c("Observed", "Population")
  }
  
  diag(adj_mat) <- 0
  for (i in 1:length(clst)) {
    x <- clst[[i]]
    adj_mat[x$area,x$area] <- adj_mat[x$area,x$area] * (10 * i)
  }
  
  setting$model <- model
  setting$stattype <- stattype
  setting$scanmethod <- scanmethod
  setting$cartesian <- !latlon
  setting$rantype <- rantype
  
  input <- list()
  input$coordinates <- coordinates
  input$case <- case
  input$adj_mat <- adj_mat
  
  retval <- list(call = call, input = input, cluster = clst,
                 setting = setting, comments = comments)
  class(retval) <- "rflexscan"

  return(retval)
}

#' Print rflexscan object
#' 
#' Print method for the rflexscan object.
#' 
#' @param x
#' An rflexscan object to be printed.
#' 
#' @param ...
#' Ignored.
#' 
#' @seealso \link{rflexscan}
#' 
#' @method print rflexscan
#' @export
#' 
print.rflexscan <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Most likely cluster (P-value: ", x$cluster[[1]]$pval, "):\n", sep = "")
  cat(x$cluster[[1]]$name, fill = 76)
  cat("Number of secondary clusters:", length(x$cluster) - 1, "\n\n")
}


#' Print rflexscanCluster object
#' 
#' Print method for the rflexscanCluster object.
#' 
#' @param x
#' An rflexscanCluster object to be printed.
#' 
#' @param ...
#' Ignored.
#' 
#' @method print rflexscanCluster
#' @export
#' 
print.rflexscanCluster <- function(x, ...) {
  cat("\n")
  cat("Areas included ...........:\n")
  cat(x$name, fill = 76)
  cat("Maximum distance .........: ", x$max_dist, "\n")
  cat("(areas: ", x$from, " to ", x$to, ")\n", sep = "")
  cat("Number of cases ..........:", x$n_case, "\n")
  
  if (!is.null(x$expected)) {
    cat("Expected number of cases .:", x$expected, "\n")
    cat("Overall relative risk ....:", x$RR, "\n")
  } else if (!is.null(x$population)) {
    cat("Population ...............:", x$population, "\n")
  }
  
  cat("Statistic value ..........:", x$stats, "\n")
  cat("Monte Carlo rank .........:", x$rank, "\n")
  cat("P-value ..................:", x$pval, "\n")
  cat("\n")
}

#' Summarizing rflexscan results
#' 
#' Summary method for rflexscan objects.
#' 
#' @param object
#' An rflexscan object to be summarized.
#' 
#' @param ...
#' Ignored.
#' 
#' @seealso \link{rflexscan}
#' 
#' @method summary rflexscan
#' @export
#' 
summary.rflexscan <- function(object, ...) {
  n_cluster <- length(object$cluster)
  total_areas <- nrow(object$input$case)
  total_cases <- sum(object$input$case[,"Observed"])

  n_area <- sapply(object$cluster, function(i){length(i$area)})
  max_dist <- sapply(object$cluster, function(i) {i$max_dist})
  n_case <- sapply(object$cluster, function(i) {i$n_case})
  stats <- sapply(object$cluster, function(i) {i$stats})
  pval <- sapply(object$cluster, function(i) {i$pval})
  
  if (toupper(object$setting$model) == "POISSON") {
    expected <- sapply(object$cluster, function(i) {i$expected})
    RR <- sapply(object$cluster, function(i) {i$RR})

    table <- data.frame(NumArea=n_area, MaxDist=max_dist, Case=n_case, 
                        Expected=expected, RR=RR, Stats=stats, P=pval)
  } else {
    population <- sapply(object$cluster, function(i) {i$population})

    table <- data.frame(NumArea=n_area, MaxDist=max_dist, Case=n_case,
                        Population=population, Stats=stats, P=pval)
  }
  row.names(table) <- 1:n_cluster

  retval <- list(call=object$call,
                 total_areas=total_areas, total_cases=total_cases,
                 cluster=table, setting=object$setting)
  
  class(retval) <- "summary.rflexscan"
  return(retval)
}



#' Print summary of flexscan results
#' 
#' Print summary of flexscan results to the terminal.
#' 
#' @param x
#' An summary.rflexscan object to be printed.
#' 
#' @param ...
#' Ignored.
#' 
#' @seealso \link{rflexscan}, \link{summary.rflexscan}
#' 
#' @method print summary.rflexscan
#' @export
#' 
print.summary.rflexscan <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Clusters:\n")
  signif <- symnum(x$cluster[,"P"], corr = FALSE, 
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                   symbols = c("***", "**", "*", ".", " "))
  
  dig <- ceiling(log10(x$setting$simcount))
  
  if (toupper(x$setting$model) == "POISSON") {
    table <- data.frame(NumArea = x$cluster$NumArea,
                   MaxDist = round(x$cluster$MaxDist, 3),
                   Case = x$cluster$Case,
                   Expected = round(x$cluster$Expected, 3),
                   RR = round(x$cluster$RR, 3),
                   Stats = round(x$cluster$Stats, 3),
                   P = format(round(x$cluster$P, dig), nsmall = dig),
                   signif)
  } else {
    table <- data.frame(NumArea = x$cluster$NumArea,
                   MaxDist = round(x$cluster$MaxDist, 3),
                   Case = x$cluster$Case,
                   Population = x$cluster$Population,
                   Stats = round(x$cluster$Stats, 3),
                   P = format(round(x$cluster$P, dig), nsmall = dig),
                   signif)
  }
  colnames(table)[ncol(table)] <- ""
  print(table, quote = FALSE, right = TRUE, print.gap = 2)
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
  
  cat("Limit length of cluster:", x$setting$clustersize, "\n")
  cat("Number of areas:", x$total_areas, "\n")
  cat("Total cases:", x$total_cases, "\n")
  if (x$setting$cartesian) {
    cat("Coordinates: Cartesian\n")
  } else {
    cat("Coordinates: Latitude/Longitude\n")
  }
  cat("Model:", x$setting$model, "\n")
  cat("Scanning method:", x$setting$scanmethod, "\n")
  cat("Statistic type:", x$setting$stattype, "\n\n")
}
  

#' Graph plotting of flexscan results
#' 
#' Display detected clusters by a graph representation.
#' 
#' @param x
#' An rflexscan object.
#' 
#' @param rank
#' An integer vector which specifies ranks of clusters to be displayed.
#' 
#' @param pval
#' A threshold of P-value. Clusters with P-values of <\code{pval} will be displayed.
#' 
#' @param vertexsize
#' Size of vertex of the graph.
#' 
#' @param xlab
#' A label of the x axis.
#' 
#' @param ylab
#' A label of the y axis.
#' 
#' @param xlim
#' The x limits of the plot.
#' 
#' @param ylim
#' The y limits of the plot.
#' 
#' @param col
#' A vector of colors for each cluster.
#' 
#' @param frame_color
#' Color of frames in the graph.
#' 
#' @param vertex_color
#' Fill color of vertices that are not included in any clusters.
#' 
#' @param ...
#' Other parameters to be passed to \link{plot.igraph} function.
#' 
#' @details 
#' Clusters are colored using the current palette. Please use \link{palette}
#' function to specify colors of each cluster. Note that clusters with ranks
#' larger than the number of colors in the palette are not highlighted.
#' 
#' @seealso \link{rflexscan}
#' 
#' @examples
#' # load sample data (North Carolina SIDS data)
#' library(spdep)
#' data("nc.sids")
#' 
#' # calculate the expected numbers of cases
#' expected <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)
#' 
#' # run FleXScan
#' fls <- rflexscan(x = nc.sids$x, y = nc.sids$y,
#'                  observed = nc.sids$SID74,
#'                  expected = expected,
#'                  name = rownames(nc.sids),
#'                  clustersize = 10,
#'                  nb = ncCR85.nb)
#' 
#' # display all clusters
#' plot(fls)
#' 
#' # display clusters with rank 1, 2 and 3
#' plot(fls, rank = c(1, 2, 3))
#' 
#' # display clusters of P-value <= 0.05
#' plot(fls, pval = 0.05)
#' 
#' @importFrom igraph graph_from_adjacency_matrix V V<- E E<- plot.igraph
#' 
#' @method plot rflexscan
#' @export
#' 
plot.rflexscan <- function(x,
                  rank=1:length(x$cluster),
                  pval=1,
                  vertexsize=max(x$input$coordinates[,1])-min(x$input$coordinates[,1]),
                  xlab=colnames(x$input$coordinates)[1],
                  ylab=colnames(x$input$coordinates)[2],
                  xlim=c(min(x$input$coordinates[,1]), max(x$input$coordinates[,1])),
                  ylim=c(min(x$input$coordinates[,2]), max(x$input$coordinates[,2])),
                  col=palette(),
                  frame_color="gray40",
                  vertex_color="white",
                  ...) {
  
  g <- graph_from_adjacency_matrix(x$input$adj_mat, mode = "undirected", diag = FALSE, weighted = TRUE)
  V(g)$size <- vertexsize
  V(g)$frame.color <- frame_color
  V(g)$color <- vertex_color
  V(g)$label <- ""
  E(g)$color <- frame_color
  
  # color clusters
  for (i in 1:min(length(col), length(x$cluster))) {
    if (i %in% rank & x$cluster[[i]]$pval <= pval) {
      V(g)$color[x$cluster[[i]]$area] <- col[i]
      E(g)$color[E(g)$weight == 10 * i] <- col[i]
    }
  }
  
  if (x$setting$cartesian) {
    plot(g, axes = TRUE, layout = as.matrix(x$input$coordinates[,c(1,2)]), rescale = FALSE,
         xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
  } else {
    # flip X-Y (x:longitude, y:latitude)
    plot(g, axes = TRUE, layout = as.matrix(x$input$coordinates[,c(2,1)]), rescale = FALSE,
         xlab = ylab, ylab = xlab, xlim = ylim, ylim = xlim, ...)
  }
}


#' Display choropleth map
#' 
#' Display choropleth map of detected clusters.
#' 
#' @param polygons
#' A SpatialPolygonsDataFrame.
#' 
#' @param fls
#' An rflexscan object.
#' 
#' @param col
#' A vector of colors for each cluster.
#' 
#' @param region_color
#' Color of regions that are not included in any clusters. 
#' 
#' @param rank
#' An integer vector which specifies ranks of clusters to be displayed.
#' 
#' @param pval
#' A threshold of P-value. Clusters with P-values of <\code{pval} will be displayed.
#' 
#' @param ...
#' Other parameters to be passed to plot function.
#' 
#' @details 
#' Clusters are colored using the current palette. Please use \link{palette}
#' function to specify colors of each cluster. Note that clusters with ranks
#' larger than the number of colors in the palette are not highlighted.
#'
#' @seealso \link{rflexscan}
#' 
#' @examples
#' \donttest{
#' # load sample data (North Carolina SIDS data)
#' library(rgdal)
#' library(spdep)
#' data("nc.sids")
#' sids.shp <- readOGR(system.file("shapes/sids.shp", package="spData")[1])
#' 
#' # calculate the expected numbers of cases
#' expected <- nc.sids$BIR74 * sum(nc.sids$SID74) / sum(nc.sids$BIR74)
#' 
#' # run FleXScan
#' fls <- rflexscan(x = nc.sids$x, y = nc.sids$y,
#'                  observed = nc.sids$SID74,
#'                  expected = expected,
#'                  name = rownames(nc.sids),
#'                  clustersize = 10,
#'                  nb = ncCR85.nb)
#' 
#' # display all clusters
#' choropleth(sids.shp, fls)
#' 
#' # display clusters with rank 1, 2 and 3
#' choropleth(sids.shp, fls, rank = c(1, 2, 3))
#' 
#' # display clusters of P-value <= 0.05
#' choropleth(sids.shp, fls, pval = 0.05)
#' }
#' 
#' @import sp grDevices graphics stats utils
#' 
#' @export
#' 
choropleth <- function(polygons,
                       fls,
                       col=palette(),
                       region_color="#F0F0F0",
                       rank=1:length(fls$cluster),
                       pval=1,
                       ...) {
  # color clusters
  for (i in 1:min(length(col), length(fls$cluster))) {
    if (!(i %in% rank & fls$cluster[[i]]$pval <= pval)) {
      col[i] <- region_color
    }
  }
  col <- c(col, region_color)

  index <- numeric(nrow(fls$input$case))
  for (i in 1:length(fls$cluster)) {
    index[fls$cluster[[i]]$area] <- i
  }
  index[index == 0 | index > length(col)] <- length(col)
  plot(polygons, col = col[index], lwd = 0.1, ...)
  box()
}
