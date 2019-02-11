#' Fit models relating sleep, age and hippocampus volume.
#'
#' @param data A dataframe.
#'
#' @return List of all fitted models.
#' @export
#'
#' @examples
#' #library(metasleep)
#' 
fit_models <- function(data){
  
  message("Checking dataset")
  data <- check_data(data)

  message("Dataset ok.")
  
  message("\n\nFitting main effect of age on sleep\n")
  
  pb <- progress_estimated(length(psqi_variables))
  age_models <- map(psqi_variables, function(y){
    pb$tick()$print()
    map(spline_bases$age_models, function(basis) {
      df <- create_design_matrix(data, list(basis), y)
      # Full model
      form <- reformulate(c("SexFemale", smoothvars(df)), response = y, intercept = TRUE)
      res <- fit_mixed_model(form, df)
      
      res[["basis"]] = list(basis)
      return(res)
    })
  })
  
  models <- map(c("Hippocampus", "TBV"), function(response){
    
    message(paste("\n\nFitting effect of sleep and age on", response, "\n"))
    pb <- progress_estimated(length(psqi_variables))
    sleep_hippocampus_models <- map(psqi_variables, function(y){
      pb$tick()$print()
      map(spline_bases$sleep_hippocampus_models[[y]], function(basis) {
        . = NULL
        df <- data %>%
          filter_at(vars(matches(paste0("^", y, "$"))), all_vars(is.na(.) == FALSE)) %>%
          create_design_matrix(basis, "")
  
        # Full model
        form <- reformulate(c("SexFemale", "ICV", smoothvars(df)), response = response,
                            intercept = TRUE)
  
        res <- c(fit_mixed_model(form, df))
  
        res[["basis"]] = list(basis)
        return(res)
      })
    })
    
  
    message(paste("\n\nFitting effect of sleep and timepoint on", response, "\n"))
    pb <- progress_estimated(length(psqi_variables))
    sleep_hippocampus_change <- map(psqi_variables, function(y){
      pb$tick()$print()
      map(spline_bases$sleep_hippocampus_change[[y]], function(basis) {
        . = NULL
        df <- data %>%
          filter_at(vars(matches(paste0("^", y, "$"))), all_vars(is.na(.) == FALSE)) %>%
          create_design_matrix(basis, "")
  
        # Full model
        form <- reformulate(c("SexFemale", "ICV", smoothvars(df)), response = response,
                            intercept = TRUE)
  
        res <- c(fit_mixed_model(form, df))
  
        res[["basis"]] = list(basis)
        return(res)
      })
    })

    message(paste("\n\nFitting interaction effect of sleep and age on", response, "\n"))
    pb <- progress_estimated(length(psqi_variables))
    sleep_hippocampus_interaction_models <- map(psqi_variables, function(y){
      pb$tick()$print()
      map(spline_bases$sleep_hippocampus_interaction_models[[y]], function(basis){
        . = NULL
        df <- data %>%
          filter_at(vars(matches(paste0("^", y, "$"))), all_vars(is.na(.) == FALSE)) %>%
          create_design_matrix(basis, "")
  
        form <- reformulate(c("SexFemale", "ICV", smoothvars(df)), response = response,
                            intercept = TRUE)
  
        res <- c(fit_mixed_model(form, df))
  
        res[["basis"]] = list(basis)
  
        return(res)
      })
    })
    
  
    message(paste("\n\nFitting interaction effect of sleep, timepoint and age on", response, "\n"))
    pb <- progress_estimated(length(psqi_variables))
    sleep_hippocampus_interaction_change <- map(psqi_variables, function(y){
      pb$tick()$print()
      map(spline_bases$sleep_hippocampus_interaction_change[[y]], function(basis){
        . = NULL
        df <- data %>%
          filter_at(vars(matches(paste0("^", y, "$"))), all_vars(is.na(.) == FALSE)) %>%
          create_design_matrix(basis, "")
  
        form <- reformulate(c("SexFemale", "ICV", smoothvars(df)), response = response,
                            intercept = TRUE)
  
        res <- c(fit_mixed_model(form, df))
  
        res[["basis"]] = list(basis)
  
        return(res)
      })
    })
  })

  list(age_models = age_models,
       sleep_hippocampus_models = models[[1]]$sleep_hippocampus_models,
       sleep_tbv_models = models[[2]]$sleep_hippocampus_models,
       sleep_hippocampus_change = models[[1]]$sleep_hippocampus_change,
       sleep_tbv_change = models[[2]]$sleep_hippocampus_change,
       sleep_hippocampus_interaction_models = models[[1]]$sleep_hippocampus_interaction_models,
       sleep_tbv_interaction_models = models[[2]]$sleep_hippocampus_interaction_models,
       sleep_hippocampus_interaction_change = models[[1]]$sleep_hippocampus_interaction_change,
       sleep_tbv_interaction_change = models[[2]]$sleep_hippocampus_interaction_change
       )
  
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

