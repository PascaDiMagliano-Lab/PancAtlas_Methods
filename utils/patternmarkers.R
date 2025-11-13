## author: Atul Deshpande
## email: adeshpande@jhu.edu
newPatternMarkers <- function(object, data = NULL, fdrThresh = 5e-2, cutoff = 0.2, sanityCheck = c("correlation","ratio")) {
  patternList <- colnames(object@featureLoadings)
  N <- length(patternList)
  geneList <- rownames(object@featureLoadings)
  
  if (is.null(fdrThresh)) fdrThresh <- 1/length(geneList)
  
  ## get total gene expression in each pattern
  totalGeneExprFromPattern <- object@featureLoadings %*% diag(colSums(object@sampleFactors))
  rownames(totalGeneExprFromPattern) <- geneList
  colnames(totalGeneExprFromPattern) <- patternList
  
  ## total gene expression in the cogaps approximate A*P
  totalGeneExprfromCoGAPS <- apply(totalGeneExprFromPattern,1,sum)
  markerGenes <- geneList
  if (!is.null(data)){
    data <- as.matrix(data)
    if (sanityCheck=="ratio"){
      totalGeneExprData <- apply(data[markerGenes,],1,sum)  
      niceGenes <- names(which(abs(totalGeneExprfromCoGAPS/totalGeneExprData-1) < cutoff))
    }  else if (sanityCheck=="correlation"){
      cogapsMat <- object@featureLoadings %*% t(object@sampleFactors)  
      exprCorrelation <- sapply(rownames(cogapsMat),function(i) cor(cogapsMat[i,],data[i,]))
      niceGenes <- names(which(exprCorrelation > cutoff))
    }
    markerGenes <- niceGenes
  }
  
  
  ## calculate marker scores as fraction of gene expression in each pattern
  markerScores <- totalGeneExprFromPattern[markerGenes,]/totalGeneExprfromCoGAPS[markerGenes]
  
  ## use mixtools to group genes into three groups - low, mid and high expression
  lambda = c(0.3,0.6,0.1)
  mu = c(1,50,500)
  sigma = c(10,50,500)
  exprLevels = tryCatch({
    mixGeneExpr <- mixtools::normalmixEM(totalGeneExprfromCoGAPS[markerGenes],lambda = lambda,mu = mu, sigma = sigma)
    exprLevels<- apply(mixGeneExpr$posterior,1,which.max)
  }, error = function(e) {
    exprLevels = tryCatch({
      mixGeneExpr <- mixtools::normalmixEM(totalGeneExprfromCoGAPS[markerGenes],lambda = c(0.7,0.3),mu = mu[2:3], sigma = sigma[2:3])
      exprLevels<- apply(mixGeneExpr$posterior,1,which.max) + 1
      return(exprLevels)
    }, error = function(e){
      return(rep(3, length(markerGenes)))
    })
  })
  
  names(exprLevels) <- markerGenes
  
  ## find mixtures of three distributions of the fractional expression for each expression level
  ## low expression
  mixScores <- means <- stds <- c()
  
  for (ii in 1:3){
    if (ii %in% unique(exprLevels))
    {
      temp <- mixtools::normalmixEM(markerScores[which(exprLevels==ii),],lambda = c(0.6,0.3,0.1),mean.constr = c(0,1/length(patternList),NA), sigma = c(.1/N,.5/N,.3))
      means <- c(means,temp$mu[2])
      stds <- c(stds,temp$sigma[2])
      temp <- matrix(temp$posterior[,3], nrow = sum(exprLevels==ii))
      rownames(temp) <- names(which(exprLevels==ii))    
      mixScores <- rbind(mixScores,temp)
    }
    else
    {
      means <- c(means,0)
      stds <- c(stds,1)
    }
  }
  
  ## calculate outlier threshold for each expression level
  #means <- c(mixScores[[1]]$mu[2],mixScores[[2]]$mu[2],mixScores[[3]]$mu[2])
  #means <- means[exprLevels]
  #stds <- c(mixScores[[1]]$sigma[2],mixScores[[2]]$sigma[2],mixScores[[3]]$sigma[2])
  #stds <- stds[exprLevels]
  score_minus_means <- markerScores - means
  zscores <- score_minus_means/stds
  pVals <- 1-pnorm(zscores)
  qVals <- qvalue::qvalue(pVals, pi0 = 1)
  qVals <- qVals$qvalues
  
  colnames(mixScores) <- patternList
  
  ## assign markers to patterns
  #markersByPattern <- lapply(patternList,function(pattern) names(which(markerScores[,pattern]>markerThresholds)))
  markersByPattern <- lapply(patternList,function(pattern) names(which(qVals[,pattern]<fdrThresh)))
  names(markersByPattern) <- patternList    
  
  ## order by markerScores
  #markersByPattern <- lapply(names(markersByPattern),function(pattern) markersByPattern[[pattern]][order(markerScores[markersByPattern[[pattern]],pattern],decreasing = T)])
  ## order by qValues
  markersByPattern <- lapply(names(markersByPattern),function(pattern) markersByPattern[[pattern]][order(qVals[markersByPattern[[pattern]],pattern])])
  names(markersByPattern) <- patternList    
  
  return(list(
    "PatternMarkers"=markersByPattern,
    "PatternMarkerScores"=markerScores,
    "PatternMixScores"=mixScores,
    "markerExprLevels"=exprLevels
  ))
}
