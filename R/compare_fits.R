#' Compare individual fits to meta analysis fit.
#'
#' @param meta_models Dataframe returned form \code{\link{combine_fits}}.
#'
#' @return Dataframe.
#' @export
#'
compare_fits <- function(meta_models){
  meta_models %>% 
    pmap_dfr(function(...){
      args <- list(...)
      
      grid <- set_up_grid(args$basis, args$Model, args$PSQI)
      
      betas <- c(args$study_betas, list(t(args$metafit$coefficients)))
      names(betas) <- c(paste0("Study", seq(1, length(args$study_betas), by = 1)), "Metamodel")
    
      map_df(betas, ~ grid$X %*% .x) %>% 
        bind_cols(grid$df) %>% 
        mutate(
          Model = args$Model,
          PSQI = args$PSQI,
          Criterion = args$Criterion
        ) %>% 
        gather(key = "Dataset", value = "Fit", starts_with("Study"), .data$Metamodel) %>% 
        group_by(.data$Model, .data$PSQI, .data$Criterion) %>% 
        nest()
    })  
}
