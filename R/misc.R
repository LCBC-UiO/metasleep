
predict_mat <- function(bs, df){
  X <- PredictMat(bs, data = df)
  # Detect the intercept column and drop it if it exists
  keep_inds <- setdiff(seq(1, ncol(X)), which(apply(X, 2, function(x) length(unique(x)) == 1)))
  X <- X[, keep_inds, drop = FALSE]
  colnames(X) <- paste0(paste(bs$term, collapse = "_"), "_bs", seq(ncol(X)))
  as_tibble(X)
}
