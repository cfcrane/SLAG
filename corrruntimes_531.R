x <- c(407.29, 80.62, 164.71, 90.81, 196.10, 104.33, 71.62, 217.33, 178.95)
y <- c(1458.00, 155.00, 509.00, 669.00, 638.00, 143.44, 220.00, 3820.62, 667.00)
capture.output(cor.test(x, y, alternative="two.sided", method="pearson"), file = "hemiStanleyatramtoslagruntimecorr.txt")
