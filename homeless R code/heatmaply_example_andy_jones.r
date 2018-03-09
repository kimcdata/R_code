







#With dimension argument, times 100 for tiff

plotHeatMapDims <- function(matrix, fileName, plotTitle, xDim, yDim){
  
  #quantile.range <- quantile(matrix, probs = seq(0, 1, 0.01))
  
  #palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  
  #color.palette  <- colorRampPalette(c("#FC8D59", "#FFFFBF", "#91CF60"))(length(palette.breaks) - 1)
  
  pdf(file=paste(fileName,".pdf",sep =""),width = xDim, height = yDim)
  
  dist.pear <- function(x) as.dist(1-cor(t(x)))
  
  hclust.ave <- function(x) hclust(x, method="average")
  
  
  
  heatmap.2(
    
    matrix,
    
    scale      = "row",
    
    trace      = "none",
    
    key        = TRUE,
    
    col    = bluered(50),
    
    breaks = seq(-2.5, 2.5, length.out = 51),
    
    main = plotTitle,
    
    margin=c(10, 10),
    
    distfun=dist.pear,
    
    hclustfun=hclust.ave
    
  )
  
  dev.off()
  
  
  
  tiffxDim = xDim * 100
  
  tiffyDim = yDim * 100
  
  tiff(file=paste(fileName,".tif",sep =""),width = tiffxDim, height = tiffyDim, units = "px", pointsize = 24,
       
       compression =  "lzw")
  
  heatmap.2(
    
    matrix,
    
    scale      = "row",
    
    trace      = "none",
    
    key        = TRUE,
    
    col    = bluered(50),
    
    breaks = seq(-2.5, 2.5, length.out = 51),
    
    
    main = plotTitle,
    
    margin=c(10, 10),
    
    distfun=dist.pear,
    
    hclustfun=hclust.ave
    
  )
  
  dev.off()
  
  hm = heatmaply(fullKinaseMatrix, scale = "row", file = "heatmaply.png", distfun = dist.pear, hclustfun = hclust.ave, width = 1600, height = 1600, margins = c(250,100,NA,0), colors = c("blue","white","red"), k_col = 2, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue",high = "red",mid = "white", midpoint = 0, limits = ))
  
  return(hm)
  
}



kinaseExpressionFrame = read.csv('kinaseExpression_FANTOM5_plus_classes_plus_TNNI3.csv')

row.names(kinaseExpressionFrame) = kinaseExpressionFrame$Kinase

kinaseExpressionFrame = kinaseExpressionFrame[5:80]

fullKinaseMatrix = data.matrix(kinaseExpressionFrame)#make data into a matrix

hm = plotHeatMapDims(fullKinaseMatrix,"all_kinases_heat_map_with_TNNI3_pearcorr",
                
                "Heat map of all kinases across human tissues with TNNI3", 30, 60)




heatmaply(fullKinaseMatrix, scale = "row", file = "heatmaply.png", distfun = dist.pear, hclustfun = hclust.ave, width = 1600, height = 5600, margins = c(250,100,NA,0), k_col = 2, scale_fill_gradient_fun = ggplot2::scale_fill_gradientn(n = 100, colours = c("blue","white","red"), values = c(seq(0, 0.2, length.out=40),seq(0.2,0.35,length.out=20),seq(0.35,1,length.out=40))))