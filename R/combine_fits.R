#' Combine fitted models
#'
#' @param models Returned from \code{\link{fit_models}}.
#'
#' @return list
#' @export
#' 
#'
combine_fits <- function(models){
  # w <- models[[1]]
  # x <- w[[2]]
  # y <- x[[1]]
  # z <- y[[1]]
  model_df <- map_dfr(models, function(w){
    # Model level
    map_df(w, function(x){
      # PSQI level
      map_df(x, function(y){
        # Knot level
        map_df(y, function(z){
          tibble(
            beta = list(z$beta),
            S = list(z$S),
            logLik = as.numeric(z$logLik),
            df = attr(z$logLik, "df"),
            nobs = attr(z$logLik, "nobs"),
            basis = z$basis,
            knots = count_knots(z$basis)
            )
          })
      }, .id = "PSQI")
    }, .id = "Model") 
  }, .id = "Dataset") %>% 
    group_by(.data$Model, .data$PSQI, .data$df, .data$knots) %>% 
    summarise(
      beta = list(.data$beta),
      S = list(.data$S),
      logLik = sum(.data$logLik),
      basis = .data$basis[1],
      nobs = sum(.data$nobs)
    ) %>% 
    ungroup() %>% 
    mutate(
      AIC = -2 * .data$logLik + 2 * .data$df,
      BIC = -2 * .data$logLik + log(.data$nobs) * .data$df
    ) %>% 
    gather(key = "Criterion", value = "value", .data$AIC, .data$BIC)
    
  # Pick the model with lowest AIC/BIC for each combination
  meta_models <- model_df %>% 
    group_by(.data$Model, .data$PSQI, .data$Criterion) %>% 
    filter(.data$value == min(.data$value)) %>% 
    ungroup() %>% 
    pmap_dfr(function(...) {
      df <- list(...) 
      study_betas <- t(do.call(cbind, df$beta)) 
      metafit <- mvmeta(study_betas ~ 1, S = df$S, method = "mm")
      as_tibble(df[c("Model", "PSQI", "knots", "Criterion", "value")]) %>%
        distinct() %>% 
        mutate(
          metafit = list(metafit), 
          basis = list(df$basis), 
          study_betas = list(df$beta),
          pval = list(df$pval))
    })
  
  
  list(meta_models = meta_models, 
       model_df = model_df)
}



set_up_grid <- function(basis, Model, PSQI){
  if(Model == "age_models"){
    df <- age_grid
    # Intercept, SexFemale
    X <- as.matrix(cbind(0, 1, predict_mat(basis, df)))
  } else if(Model == "sleep_hippocampus_models" || Model == "sleep_hippocampus_interaction_models") {
    df <- crossing(age_grid, psqi_grid[[PSQI]])
    # Intercept, SexFemale, ICV
    X <- cbind(0, 1, 0, as.matrix(map_dfc(basis, ~ predict_mat(.x, df))))
    df <- rename(df, PSQI_value = !!sym(PSQI))
  } else {
    stop("Unknown model.")
  }
  list(X = X, df = df)
}


#' Compute confidence bands for the fits
#'
#' @param meta_models Dataframe.
#'
#' @return Dataframe.
#' @export
#'
compute_confidence_bands <- function(meta_models){
  meta_models %>%
    pmap_dfr(function(...){
      args <- list(...)
      grid <- set_up_grid(args$basis, args$Model, args$PSQI) 
      betas <- rmvnorm(1000, mean = c(args$metafit$coefficients), as.matrix(args$metafit$vcov))
      confints <- t(apply(grid$X %*% t(betas), 1, function(x) c(mean(x), quantile(x, probs = c(0.025, 0.975), names = FALSE))))
      colnames(confints) <- c("fit", "lower", "upper")
      as_tibble(confints) %>% 
        bind_cols(grid$df) %>% 
        mutate(
          Model = args$Model,
          PSQI = args$PSQI,
          Criterion = args$Criterion
        )
    })
}


count_knots <- function(bs){
    if(inherits(bs, "mgcv.smooth")){
      if("margin" %in% names(bs)){
        res <- paste0("(", paste(map(bs$margin, ~ .x$bs.dim), collapse = ", "), ")")
      } else {
        res <- bs$bs.dim
      }
    } else {
      res <- map(bs, count_knots)
    }
  paste(res, collapse = ", ")
  }
  