df <- read.table("primermeansforR.txt")
H <- df$percent_working
N <- df$coverage
#H <- c(91, 86.8, 77.1, 54, 91.4, 82.6, 67.2, 77.6)
#N <- c("63.4x", "31.6x", "15.8x", "7.9x", "63.4x", "31.6x", "15.8x", "7.9x")
pdf(file="basicbarplotofworkingprimers.pdf", height = 4, width = 4)
barplot(H, names.arg=N, xlab="aTRAM                       SLAG", ylab="percent working primers", col="#555555", main="Primer Result", border="#555555", cex.names=0.6)
dev.off()
