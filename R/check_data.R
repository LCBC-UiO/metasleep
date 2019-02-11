check_data <- function(data){
  
  coldiff <- setdiff(names(data_template), names(data))
  if(length(coldiff) > 0){
    stop(paste("Missing variable(s) in dataframe:", paste(coldiff, collapse = ", ")))  
  }
  
  # Check that there are no missing IDs
  if(data %>% filter(is.na(.data$ID)) %>% nrow() > 0){
    stop("There are IDs with value NA. Cannot continue.")
  }
  
  # Checking the type of ID  
  if(!is.factor(pull(data, .data$ID))){
    if(is.numeric(pull(data, .data$ID)) || is.character(pull(data, .data$ID))){
      data <- data %>% 
        mutate(ID = as.factor(as.character(.data$ID)))
    } else {
      stop("ID variable must be factor, numeric, or character.")
    }
  }
  
  # Checking that IDs are not duplicated
  check <- data %>% 
    group_by(.data$ID, .data$Age) %>% 
    summarise(counts = n()) %>% 
    filter(.data$counts > 1)

  if(nrow(check) > 0){
    stop(paste("There are duplicate rows for given ID, Age pair. For IDs", 
               paste(pull(check, .data$ID), collapse = ", ")))
  }
  
  if(! "Age_bl" %in% names(data)){
    message("Adding Age_bl and Timepoint variables.")
    data <- data %>% 
      group_by(.data$ID) %>% 
      mutate(
        Age_bl = min(.data$Age),
        Timepoint = .data$Age - .data$Age_bl
        ) %>% 
      ungroup()
  }
  
  if(! "PSQI_Global" %in% names(data)){
    message("Adding PSQI_Global variable.")
    
    PSQI_mat <- data %>% 
      select_at(vars(matches("^PSQI\\_Comp[1-7]"))) %>% 
      as.matrix()
    
    PSQI_sum <- bind_cols(psqi_sum = rowSums(PSQI_mat, na.rm = TRUE), num_missing = rowSums(is.na(PSQI_mat))) %>% 
      mutate(psqi_sum = if_else(.data$num_missing <= 2, .data$psqi_sum * 7 / (7 - .data$num_missing), NA_real_)) %>% 
      pull(.data$psqi_sum)
    
    data <- data %>% 
      mutate(PSQI_Global = PSQI_sum)
  }
  
  
  message("Scaling and shifting the ICV variable for numerical stability.")
  # Using rough values from LCBC. This is important, since all datasets must use the same scaling, and not their own scalings.
  data <- data %>% 
    mutate(ICV = (.data$ICV - 1.5e6)/1.7e5)
  
  
  return(data)
  
  
}
