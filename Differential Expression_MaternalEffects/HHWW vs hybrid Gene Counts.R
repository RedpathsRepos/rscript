# top numbers are from consolidated

HHWW.down  <- c(3147, 2248, 70, 191, 91)    # matrix order equal k,m,o,p,v
HHWW.up    <- c(2646, 2653, 60, 195, 88)

hybrid.down <- c(430, 1060, 28, 83, 171)    # individual 19 from matrix O excluded - was a crazy outlier - look into more?
hybrid.up   <- c(490, 797, 50, 65, 179)


hw = (c(HHWW.down, HHWW.up))
hybrid = (c(hybrid.down, hybrid.up))


HHWW.down  <- c(9891, 8380, 253, 236, 88)    # matrix order equal k,m,o,p,v
HHWW.up    <- c(9315, 6582, 229, 645, 89)

hybrid.down <- c(1358, 3420, 159, 61, 260)    # individual 19 from matrix O excluded - was a crazy outlier - look into more?
hybrid.up   <- c(1443, 2273, 136, 110, 129)


hw = HHWW.down + HHWW.up
hybrid = hybrid.down + hybrid.up

barplot(c(mean(hw), mean(hybrid)), width = 0.1, space = .75, names.arg = c("HH vs. WW", "HW vs. WH"), ylab = "Average Number of DE Genes", xlab = "Comparison", col = c(rgb(55,126,184, maxColorValue = 255), rgb(228,26,28, maxColorValue = 255)))
abline(0,0)

chisq.test(hw, hybrid)



