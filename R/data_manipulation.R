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

## phangorn::phyDat has inexplicably stopped working.  This function 'fixes' it
## For simplicity I have not retained support for contrast matrices or ambiguity.
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

Min <- function (x, inapp.level) {
  if (length(inapp.level)) return(sum(2^(c(0:(inapp.level-2), inapp.level:12)) %in% unique(x)))
  return (sum(2^(0:12) %in% unique(x)))
}

PrepareDataFitch <- function (data) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(data)
  nam <- at$names
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(data, "contrast")
  nTip <- length(data)
  info.amounts <- InfoAmounts(data)
  max.length <- max(vapply(info.amounts, length, integer(1)))
  info <- array(0, dim=c(nChar, max.length))
  for (i in 1:nChar) info[i, 1:length(info.amounts[[i]])] <- info.amounts[[i]]
  at$names <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  data <- unlist(data, FALSE, FALSE)
  ret <- tmp[data] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inapp.level <- which(at$levels == "-")
  attr(ret, 'inapp.level') <- 2 ^ (inapp.level - 1)
  attr(ret, 'dim') <- c(nChar, nTip)  
  attr(ret, 'unique.tokens') <- apply(ret, 1, function(x) Min(x, inapp.level))
  applicable.tokens <- setdiff(powers.of.2, 2 ^ (inapp.level - 1))
  attr(ret, 'split.sizes') <- apply(ret, 1, function(x) vapply(applicable.tokens, function (y) sum(x == y), integer(1)))
  attr(ret, 'info.amounts') <- info
  dimnames(ret) <- list(NULL, nam)
  class(ret) <- 'fitchDat'
  ret
}