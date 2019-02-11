library(tidyverse)
library(mgcv)
# Name of sleep variables.
psqi_variables <- c("PSQI_Comp1_SleepQuality", "PSQI_Comp2_Latency", "PSQI_Comp3_Duration", 
  "PSQI_Comp4_Efficiency", "PSQI_Comp5_Problems", "PSQI_Comp6_Medication", 
  "PSQI_Comp7_Tired", "PSQI_Global")
names(psqi_variables) <- psqi_variables

psqi_grid <- map_dfr(psqi_variables, 
                     ~ tibble(key = .x, 
                              value = if_else(.x == "PSQI_Global", list(seq(0, 21, by = 8)), list(seq(0, 3))),
                              num_knots = if_else(.x == "PSQI_Global", list(seq(2, 3)), list(seq(2, 3)))
                              ))

# Age grid over which to plot predictions and define knots
age_grid <- tibble(Age = seq(15, 90, length.out = 5))
age_bl_grid <- tibble(Age_bl = seq(15, 90, length.out = 5))
# Timepoint grid. Subject to change
timepoint_grid <- tibble(Timepoint = seq(0, 7, length.out = 3))

# Main effect of age on PSQI variable, s(Age)
age_models <- map(seq(2, 4, by = 1), function(k){
  smoothCon(s(Age, k = k, fx = TRUE, bs = if_else(k > 2, "tp", "re")), 
            data = age_grid)[[1]]  
})

# Setting up model with s(PSQI) + s(Age) and s(PSQI) + s(Age) + ti(Age, PSQI)
# First basis for s(PSQI)
psqi_bs <- pmap(psqi_grid, function(key, value, num_knots) {
  map(num_knots, function(k) {
    smooth <- paste0("s(", key, ", k = k, fx = TRUE, bs = ", if_else(k > 2, "'tp'", "'re'"), ")")
    smoothCon(eval(rlang::parse_expr(smooth)),  data = tibble(!!key := value))[[1]]
    })
  })

# Then, basis for s(Age)
age_bs <- map(seq(2, 4, by = 1), function(k){
  smoothCon(s(Age, k = k, fx = TRUE, bs = if_else(k > 2, "tp", "re")), 
            data = tibble(Age = seq(15, 90)))[[1]]
})

# Then the interaction basis
psqi_age_bs <- pmap(psqi_grid, function(key, value, num_knots){
  grid <- crossing(age_grid, tibble(!!key := value))
  
  cross(list(k_age = seq(2, 3, by = 1), k_psqi = seq(2, 2, by = 1))) %>% 
    map(function(x){
      bs <- paste0("c(", paste(if_else(c(x$k_age, x$k_psqi) > 2, "'tp'", "'re'"), collapse = ", "), ")")
      k <- paste0("c(", max(3, x$k_age), ", ", max(3, x$k_psqi), ")")
      smooth <- paste0("ti(Age, ", key, ", k = ", k, ", bs = ", bs, ")")
      smoothCon(eval(rlang::parse_expr(smooth)), data = grid)[[1]]
    } )
})

# Generate all combinations
sleep_hippocampus_models <- map(psqi_bs, ~ cross2(.x, age_bs))
names(sleep_hippocampus_models) <- psqi_variables

sleep_hippocampus_interaction_models <- map2(psqi_bs, psqi_age_bs, ~ cross3(.x, .y, age_bs))
names(sleep_hippocampus_interaction_models) <- psqi_variables

# Finally, the Age_bl and Timepoint models
age_bl_bs <- map(seq(2, 4, by = 1), function(k){
  smoothCon(s(Age_bl, k = k, fx = TRUE, bs = if_else(k > 2, "tp", "re")), 
            data = age_bl_grid)[[1]]
})

timepoint_bs <- map(seq(2, 2, by = 1), function(k){
  smoothCon(s(Timepoint, k = k, fx = TRUE, bs = if_else(k > 2, "tp", "re")), 
            data = timepoint_grid)[[1]]
})

age_bl_timepoint_bs <- cross(list(k_age_bl = seq(2, 3, by = 1), k_tp = seq(2, 2, by = 1))) %>% 
    map(function(x){
      bs <- paste0("c(", paste(if_else(c(x$k_age_bl, x$k_tp) > 2, "'tp'", "'re'"), collapse = ", "), ")")
      k <- paste0("c(", max(3, x$k_age_bl), ", ", max(3, x$k_tp), ")")
      smooth <- paste0("ti(Age_bl, Timepoint, k = ", k, ", bs = ", bs, ")")
      smoothCon(eval(rlang::parse_expr(smooth)), data = crossing(age_bl_grid, timepoint_grid))[[1]]
    } )

sleep_hippocampus_change <- map(psqi_bs, ~ cross(list(.x, age_bl_bs, timepoint_bs, age_bl_timepoint_bs)))
names(sleep_hippocampus_change) <- psqi_variables

# s(Age_bl) + s(Timepoint) + s(PSQI) + ti(Age_bl, PSQI) + ti(Age_bl, Timepoint) + ti(PSQI, Timepoint) + ti(PSQI, Age_bl, Timepoint)
psqi_age_bl_bs <- pmap(psqi_grid, function(key, value, num_knots){
  grid <- crossing(age_bl_grid, tibble(!!key := value))
  
  cross(list(k_age_bl = seq(2, 3, by = 1), k_psqi = seq(2, 2, by = 1))) %>% 
    map(function(x){
      bs <- paste0("c(", paste(if_else(c(x$k_age_bl, x$k_psqi) > 2, "'tp'", "'re'"), collapse = ", "), ")")
      k <- paste0("c(", max(3, x$k_age_bl), ", ", max(3, x$k_psqi), ")")
      smooth <- paste0("ti(Age_bl, ", key, ", k = ", k, ", bs = ", bs, ")")
      smoothCon(eval(rlang::parse_expr(smooth)), data = grid)[[1]]
    } )
})

psqi_timepoint_bs <- pmap(psqi_grid, function(key, value, num_knots){
  grid <- crossing(timepoint_grid, tibble(!!key := value))
  
  cross(list(k_timepoint = seq(2, 3, by = 1), k_psqi = seq(2, 2, by = 1))) %>% 
    map(function(x){
      bs <- paste0("c(", paste(if_else(c(x$k_timepoint, x$k_psqi) > 2, "'tp'", "'re'"), collapse = ", "), ")")
      k <- paste0("c(", max(3, x$k_timepoint), ", ", max(3, x$k_psqi), ")")
      smooth <- paste0("ti(Timepoint, ", key, ", k = ", k, ", bs = ", bs, ")")
      smoothCon(eval(rlang::parse_expr(smooth)), data = grid)[[1]]
    } )
})

psqi_age_bl_timepoint_bs <- pmap(psqi_grid, function(key, value, num_knots){
  grid <- crossing(timepoint_grid, age_bl_grid, tibble(!!key := value))
  
  cross(list(k_timepoint = seq(2, 2, by = 1), k_age_bl = seq(2, 3, by = 1), k_psqi = seq(2, 2, by = 1))) %>% 
    map(function(x){
      bs <- paste0("c(", paste(if_else(c(x$k_timepoint, x$k_age_bl, x$k_psqi) > 2, "'tp'", "'re'"), collapse = ", "), ")")
      k <- paste0("c(", max(3, x$k_timepoint), ", ", max(3, x$k_age_bl), ", ", max(3, x$k_psqi), ")")
      smooth <- paste0("ti(Timepoint, Age_bl, ", key, ", k = ", k, ", bs = ", bs, ")")
      smoothCon(eval(rlang::parse_expr(smooth)), data = grid)[[1]]
    } )
})

# Now we need to cross these bases
# s(Age_bl) + s(Timepoint) + s(PSQI) + ti(Age_bl, PSQI) + ti(Age_bl, Timepoint) + ti(PSQI, Timepoint) + ti(PSQI, Age_bl, Timepoint)
# Everything that involves PSQI must first be concatenated and crossed.
all_psqi_bs <- pmap(list(psqi_bs, psqi_age_bl_bs, psqi_timepoint_bs, psqi_age_bl_timepoint_bs), 
                    function(...) cross(list(...)))

sleep_hippocampus_interaction_change <- list()
for(i in seq_along(all_psqi_bs)){
  tmp <- all_psqi_bs[[i]]
  tmp <- rep(tmp, length(age_bl_bs))
  
  for(j in seq_along(tmp)){
    k <- floor((j - 1) / length(tmp) * length(age_bl_bs)) + 1
    tmp[[j]] <- c(tmp[[j]], age_bl_bs[k])  
  }
  
  tmp <- rep(tmp, length(timepoint_bs))
  
  for(j in seq_along(tmp)){
    k <- floor((j - 1) / length(tmp) * length(timepoint_bs)) + 1
    tmp[[j]] <- c(tmp[[j]], timepoint_bs[k])
  }
  
  tmp <- rep(tmp, length(age_bl_timepoint_bs))

  for(j in seq_along(tmp)){
    k <- floor((j - 1) / length(tmp) * length(age_bl_timepoint_bs)) + 1
    tmp[[j]] <- c(tmp[[j]], age_bl_timepoint_bs[k])
  }

  sleep_hippocampus_interaction_change[[i]] <- tmp
}

names(sleep_hippocampus_interaction_change) <- psqi_variables


spline_bases <- list(
  age_models = age_models,
  sleep_hippocampus_models = sleep_hippocampus_models,
  sleep_hippocampus_change = sleep_hippocampus_change,
  sleep_hippocampus_interaction_models = sleep_hippocampus_interaction_models,
  sleep_hippocampus_interaction_change = sleep_hippocampus_interaction_change
)



# Reassign psqi_grid as a list of tibbles
psqi_grid <- map_dfr(psqi_variables, 
                     ~ tibble(key = .x, 
                              value = if_else(.x == "PSQI_Global", list(seq(0, 21, by = 1)), list(seq(0, 3))),
                              num_knots = if_else(.x == "PSQI_Global", list(seq(2, 3)), list(seq(2, 3)))
                     ))
psqi_grid <- pmap(psqi_grid, ~ tibble(!!..1 := ..2))
names(psqi_grid) <- psqi_variables

# Age grid over which to plot predictions and define knots
age_grid <- tibble(Age = seq(15, 90, by = 1))
age_bl_grid <- tibble(Age_bl = seq(15, 90, by = 1))
# Timepoint grid. Subject to change
timepoint_grid <- tibble(Timepoint = seq(0, 7, by = ))



# Model dataset
data_template <- tibble(
  ID = factor(),
  Sex = character(),
  Age = numeric(),
  Hippocampus = numeric(),
  ICV = numeric(),
  PSQI_Comp1_SleepQuality = numeric(),
  PSQI_Comp2_Latency = numeric(),
  PSQI_Comp3_Duration = numeric(),
  PSQI_Comp4_Efficiency = numeric(),
  PSQI_Comp5_Problems = numeric(),
  PSQI_Comp6_Medication = numeric(),
  PSQI_Comp7_Tired = numeric()
)



usethis::use_data(psqi_variables, spline_bases, age_grid, age_bl_grid, timepoint_grid, 
                  psqi_grid, data_template, internal = TRUE, overwrite = TRUE)
