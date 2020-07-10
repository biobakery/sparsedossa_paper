cor2 <- function(x, y = x, use.zero = "both", method = "spearman", 
                 R = 30) {
  if(nrow(x) != nrow(y))
    stop("Number of rows between x and y must agree!")
  
  # If use all zeros then just regular correlation
  if(use.zero == "both") {
    return(cor(x, y, method = method))
  }
  
  ind_x <- x != 0
  ind_y <- y != 0
  # If randomize then replace zeros with randomized small values
  if(use.zero == "random") {
    minMeasure_x <- min(setdiff(abs(x), 0)) / 2
    minMeasure_y <- min(setdiff(abs(x), 0)) / 2
    x_fill <- x
    y_fill <- y
    l_cor <- lapply(1:R, function(r) {
      x_fill[!ind_x] <- runif(n = sum(!ind_x), min = -minMeasure_x, max=minMeasure_x)
      y_fill[!ind_y] <- runif(n = sum(!ind_y), min = -minMeasure_y, max=minMeasure_y)
      return(cor(x_fill, y_fill, method = method))
    })
    return(Reduce("+", l_cor)/R)
  }
  
  # If use either, or neither, subset iteratively
  sapply(1:ncol(x),
         function(i) {
           sapply(1:ncol(y),
                  function(j) {
                    if(use.zero == "either")
                      ind <- ind_x[, i] | ind_y[, j]
                    if(use.zero == "neither")
                      ind <- ind_x[, i] & ind_y[, j]
                    if(sum(ind) < 3) return(0)
                    return(cor(x[ind, i], y[ind, j], method = method))
                  })
         })
}
