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
    sleep_hippocampus_change <- tryCatch({
      map(psqi_variables, function(y){
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
    },
    error = function(e) return("Failed.")
    )

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
    sleep_hippocampus_interaction_change <- tryCatch({
      map(psqi_variables, function(y){
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
    },
    error = function(e) return("Failed.")
    )
    return(list(
      sleep_hippocampus_models = sleep_hippocampus_models,
      sleep_hippocampus_change = sleep_hippocampus_change,
      sleep_hippocampus_interaction_models = sleep_hippocampus_interaction_models,
      sleep_hippocampus_interaction_change = sleep_hippocampus_interaction_change
    ))
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




