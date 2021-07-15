rdirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  logx <- log(matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE))
  # normalize by per-sample min to avoid Inf
  logx <- logx - apply(logx, 1, min) - 745
  
  sm <- apply(logx, 1, function(x) sum(exp(x)))
  return(exp(logx) / sm)
}