counts <- cpm(y.normalized, normalized.lib.sizes=TRUE)
pc <- prcomp(counts)

#plot(pc$rotation[,1], pc$rotation[,2])
plot(pc$rotation[1:11, 1], pc$rotation[1:11, 2], pch=21, bg="orange")
points(pc$rotation[12:22, 1], pc$rotation[12:22, 2], pch=21, bg="blue")