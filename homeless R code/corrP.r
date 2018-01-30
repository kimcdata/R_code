corrP <- function(table, xlab, ylab){

plot(table[,1], table[,2], pch=19, xlab=xlab, ylab=ylab)

p <- t.test(table[,1], table[,2])$p.value
cor <- cor(table[,1], table[,2], method="spearman")

return(c(p, cor))


}

