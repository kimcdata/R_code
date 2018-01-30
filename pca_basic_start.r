########################### BASIC PRINCIPLE COMPONENT ANALYSIS IN R ###############################

library(RColorBrewer)

numberGroups <- #fill this in
sampleGroupLabels <- #fill this in
sampleLabels <- # individual sample names/timepoints/labels, fill this in
yourData <- #your matrix of data

pal <- brewer.pal(n = numberGroups, name = "Set1") #see RColorBrewer help for colour palettes/sets

prcomp(t(yourData)) -> pca
summary(pca)$importance[,"PC1"]["Proportion of Variance"]*100 -> pc1Prop
summary(pca)$importance[,"PC2"]["Proportion of Variance"]*100 -> pc2Prop

XaxisRange <- c(-50,50) # Change these to make sure everything fits in your plot
YaxisRange <- c(-50,50) # Change these to make sure everything fits in your plot

x11()

plot(pca$x, pch=19, cex=1.5, col=pal[factor(sampleGroupLabels)], xlim=c(-45,45), ylim=c(-25,45), xlab=paste("PC1: ",round(pc1Prop,2),"%",sep=""),ylab=c(paste("PC2: ",round(pc2Prop,2),"%",sep="")))

# Run this if you want sample labels on your plot
text(x = pca$x[,1]+1.8, y=pca$x[,2]+1.8, labels = sampleLabels)

# Run this if you want a legend showing your sample groups

legendPos <- "top" #where your legend will be, options "top","topright","topleft","bottomleft", etc
legend(legendPos, legend = unique(sampleGroupLabels), pch=19, col=pal[factor(unique(sampleGroupLabels))])
       
########################### EXTRACTING GENES LINKED TO A SPECIFIC PRINCIPLE COMPONENT #########################

prcomp(t(yourData)) -> pca
pca$