qcut = function(x, n) {
  quantiles = seq(0, 1, length.out = n+1)
  cutpoints = unname(quantile(x, quantiles, na.rm = TRUE))
  as.character(cut(x, cutpoints, include.lowest = TRUE, labels=FALSE))
}