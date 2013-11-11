# values from consolidated gene counts.

WWHW.dn   <- c(248, 256, 14, 36, 50)    # matrix order equal k,m,o,p,v
WWHW.up   <- c(242, 135, 30, 39, 74)

WWWH.dn   <- c(305, 81, 15, 10, 21)    # individual 19 from matrix O excluded - was a crazy outlier - look into more?
WWWH.up   <- c(292, 139, 17, 22, 38)

HHHW.dn   <- c(1438, 105, 24, 14, 61)
HHHW.up   <- c(1221, 399, 37, 28, 43)

HHWH.dn   <- c(1109, 121, 19, 133, 24)
HHWH.up   <- c(1172, 104, 21, 157, 31)


different.mother.type <- WWHW.dn + WWHW.up
same.mother.type <- WWWH.dn + WWWH.up

different.mother.type <- HHWH.dn + HHWH.up
same.mother.type <- HHHW.dn + HHHW.up

different.mother.type <- WWHW.dn + WWHW.up + HHWH.dn + HHWH.up
same.mother.type <- WWWH.dn + WWWH.up + HHHW.dn + HHHW.up



plot(1:5, different.mother.type, ylim = c(0, 6700))
points(1:5, same.mother.type, ylim = c(0, 3500), col = "green")

data = as.matrix(t(cbind(different.mother.type, same.mother.type)))

barplot(data[], beside = TRUE, ylab = "Number of DE Genes", xlab = "Matrix",names.arg = c("1","2","3","4","5"), col = c(rgb(255,127,0, maxColorValue = 255),rgb(152,78,163, maxColorValue = 255)), legend = c("Different Mother", "Same Mother"))
abline(0,0)

chisq.test(hw, hybrid)