x <- c(1312.67, 223.05, 414.52, 389.29, 978.90, 582.81, 346.67, 215.52, 1777.62, 531.05)
y <- c(252.43, 793.50, 308.25, 230.57, 460.17, 205.48, 1512.31, 1836.75, 400.89, 7966.31)
capture.output(cor.test(x, y, alternative="two.sided", method="pearson"), file = "Stanleyatramtoslagruntimecorr.txt")
