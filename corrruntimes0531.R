x <- c(3589.90, 1609.24, 1720.48, 2494.33, 1827.24, 4696.67, 1658.00, 1750.10, 2657.57, 2055.24)
y <- c(855.33, 3049.90, 2023.14, 2583.86, 344.33, 496.86, 103.11, 2292.52, 1082.00, 122.00)
capture.output(cor.test(x, y, alternative="two.sided", method="pearson"), file = "Zeaatramtoslagruntimecorr.txt")
