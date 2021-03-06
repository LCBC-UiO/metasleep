
predict_mat <- function(bs, df){
  X <- PredictMat(bs, data = df)
  # Detect the intercept column and drop it if it exists
  keep_inds <- setdiff(seq(1, ncol(X)), which(apply(X, 2, function(x) length(unique(x)) == 1)))
  X <- X[, keep_inds, drop = FALSE]
  colnames(X) <- paste0(paste(bs$term, collapse = "_"), "_bs", seq(ncol(X)))
  as_tibble(X)
}


create_design_matrix <- function(df, basis, psqi_var){
  df %>% 
    mutate(SexFemale = as.integer(.data$Sex == "Female")) %>% 
    select_at(vars(matches(paste0("^", psqi_var, "$")), .data$ID, 
                   .data$SexFemale, .data$ICV, .data$Hippocampus, .data$TBV)) %>% 
    bind_cols(map_dfc(basis, ~ predict_mat(.x, df = df))) %>% 
    na.omit()
}

smoothvars <- function(df){
  str_subset(names(df), "\\_bs[:digit:]$")  
}


fit_mixed_model <- function(form, df){
  tryCatch({
    mod <- suppressWarnings(gamm4(form, data = df, random = ~ (1|ID), REML = FALSE)$mer)  
    
    list(
      beta = coef(summary(mod))[, "Estimate", drop = FALSE],
      S = as.matrix(vcov(mod)),
      logLik = logLik(mod)
    )
  },
  error = function(e) list(beta = NULL, S = NULL, logLik = NULL)
  )
}

