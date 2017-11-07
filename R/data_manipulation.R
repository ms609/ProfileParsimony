## Copied from phangorn file phyDat.R.  Edited for style only.
FastTable <- function (data) {                                                                                 
  if(!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)                    
  }
  da <- do.call("paste", c(data, sep = "\r"))
  ind <- !duplicated(da)
  levels <- da[ind]
  cat <- factor(da, levels = levels)
  nl <- length(levels(cat))
  bin <- (as.integer(cat) - 1)                            
  pd <- nl
  bin <- bin[!is.na(bin)]
  if (length(bin)) bin <- bin + 1
  y <- tabulate(bin, pd)
  result <- list(index = bin, weights = y, data = data[ind,])
  result                                                                              
}   

#' Load phyDat object
#'
#' A convenient wrapper for \pkg{phangorn}'s \code{phyDat}
#'
#' @param data data table, perhaps from read.nexus.data
#' @param levels tokens - values that all characters migkt take
#' @param compress Compress identical transformation series into a single row of the phyDat object
#' For simplicity I have not retained support for contrast matrices or ambiguity.
#' @return a \code{phyDat} object
#' @export
PhyDat <- function (data, levels = NULL, compress = TRUE, ...) {
  if (is.null(levels)) stop("Levels not supplied")
  nam <- names(data)
  # data <- as.data.frame(t(data), stringsAsFactors = FALSE)
  if (length(data[[1]]) == 1) {
      compress <- FALSE
  }
  if (compress) {
    ddd    <- FastTable(data)
    data   <- ddd$data
    weight <- ddd$weight
    index  <- ddd$index
    n.rows <- length(data[[1]])
  } else {
    n.rows <- length(data[[1]])
    weight <- rep(1, n.rows)
    index <- 1:n.rows
  }
  n.levels <- length(levels)
  contrast <- diag(n.levels)
  all.levels <- levels
  att <- attributes(data)
  data <- lapply(data, match, all.levels)
  attributes(data) <- att
  row.names(data) <- as.character(1:n.rows)
  data <- na.omit(data)
  tmp  <- match(index, attr(data, "na.action"))
  index <- index[is.na(tmp)]
  index <- match(index, unique(index))
  rn <- as.numeric(rownames(data))
  attr(data, "na.action") <- NULL
  weight <- weight[rn]
  n.rows <- dim(data)[1]
  names(data) <- nam
  attr(data, "row.names") <- NULL
  attr(data, "weight")    <- weight
  attr(data, "nr")        <- n.rows
  attr(data, "nc")        <- length(levels)
  attr(data, "index")     <- index
  attr(data, "levels")    <- levels
  attr(data, "allLevels") <- all.levels
  attr(data, "type")      <- "USER"
  attr(data, "contrast")  <- contrast
  class(data)             <- "phyDat"
  data
}

#' Prepare data for Profile Parsimony
#' 
#' @param data dataset of class \code{phyDat}
#' @param precision number of random trees to generate when calculating Profile curves. 
#'                  With 22 tokens (taxa), a precision increase of 4e+05 to 8e+05 sees
#'                  little difference in most steps; 80% of the differences are within 
#'                  0.05 bits, and c. 95% within 0.20 bits
#'                  values > 1e+06 consume enough memory to cause my 2012 desktop to struggle
#'
#' @return a dataset of class 'profileDat'
#'
#' @author Martin R. Smith; written with reference to phangorn:::prepareDataFitch
#' @export
PrepareDataProfile <- function (data, precision = 400000) {
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  
  at$names <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  data <- unlist(data, recursive=FALSE, use.names=FALSE)
  ret <- tmp[data]
  ret <- as.integer(ret)
  attributes(ret) <- at
  inappLevel <- which(at$levels == "-")
  attr(ret, 'inappLevel') <- 2 ^ (inappLevel - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  attr(ret, 'split.sizes') <- apply(ret, 1, function(x) vapply(applicableTokens, function (y) sum(x == y), integer(1)))
  attr(ret, 'info.amounts') <- InfoAmounts(ret, precision)
  attr(ret, 'bootstrap') <- c('info.amounts', 'split.sizes')
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- 'profileDat'
  ret
}
