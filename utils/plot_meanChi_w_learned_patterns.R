plot_cogaps_metrics <- function(cogaps_results, sample_name){
  meanChi <- lapply(cogaps_results, function(x){
    x@metadata$meanChiSq
  }) %>% unlist()
  
  meanChi <- data.frame(meanChiSquare = meanChi) %>%
    rownames_to_column(var = 'Patterns_Requested') %>%
    mutate(Patterns_Requested = as.numeric(gsub("_patterns$","",Patterns_Requested)))
  
  p1 <- ggplot(meanChi, aes(x = Patterns_Requested, y = meanChiSquare)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_continuous(breaks = seq(0,20))  +
    ggtitle(sample_name)
  
  patterns_learned <- lapply(cogaps_results, function(x){
    dim(x@sampleFactors)[2]
  }) %>% unlist()
  
  patterns_learned <- data.frame(Patterns_Learned = patterns_learned) %>%
    rownames_to_column(var = 'Patterns_Requested') %>%
    mutate(Patterns_Requested = as.numeric(gsub("_patterns$","",Patterns_Requested)))
  
  p2 <- ggplot(patterns_learned, aes(x = Patterns_Requested, y = Patterns_Learned)) +
    geom_point() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
    scale_x_continuous(breaks = seq(0,20))
  
  p1 + p2
}
