## Problems
plotInt <- function(int) {
  if (length(int) > 0) {
    for (i in 1:length(int)) {
      plot(int[[i]][[1]], add=T, pch=19, col="green")
    }
  }
}

# self intersect
x <- c(1, 3, 1, 3)
y <- c(1, 3, 3, 1)
coords <- data.frame(x, y)
l <- createLine(coords)
#zerodist(createPoint(coords))
plot(l, type="b")

int <- findSelfIntersections(sl = l)
plotInt(int)


# self intersect with only one in between node (rare)
x <- c(1, 3, 2, 3)
y <- c(1, 3, 2, 1)
coords <- data.frame(x, y)
l <- createLine(coords)
plot(l, type="b")

# intersection is a line. fix me
int <- findSelfIntersections(sl = l)
plotInt(int)

# self intersect on node. (rare)
x <- c(1, 2, 3, 1, 2,  3)
y <- c(1, 2, 3, 3, 2, 1)
coords <- data.frame(x, y)
l <- createLine(coords)
#zerodist(createPoint(coords))
plot(l, type="b")

# fix me. test if their exist a double coordinate pair
int <- findSelfIntersections(sl = l)
plotInt(int)

# self intersect on node with only one in between node (rare)
x <- c(1, 2,  3, 2, 3)
y <- c(1, 2, 3, 2,  1)
coords <- data.frame(x, y)
l <- createLine(coords)

plot(l, type="b")
# intersection is line, problaly fixed when first searched for double coordinate paires
int <- findSelfIntersections(sl = l)
plotInt(int)

###############


cc <- getCoordinatesOfALine(l)
pts <- zerodist(obj = createPoint(cc))

if (length(pts) > 2) {
  # more spikes
  # use only first part
  pts <- pts[1:2]
  
}
dp <- max(pts)
cc <- cc[1:dp, ]


cc <- cc[, ]

nl <- createLine(as.data.frame(cc))
plot(nl, add=T, col="red")


if (gIntersects(l, l)) {
  g <- gIntersection(l, )l)
  plot(g, add=T, pch=19, col="green")
}

gIsValid(l, reason = T)

gCrosses(l)
gNode(l)
gRelate(l)
gTouches(l)
gOverlaps(l)
gDisjoint(l)

#Po1 <- readWKT("POLYGON((0 0,0 2,1 3.5,3 3,4 1,3 0,0 0))")
Po1 <- readWKT("LINESTRING(0 0,0 2,1 3.5,3 3,4 1,3 0)")
Lin1 <- readWKT("LINESTRING(0 3,1 1,2 2,3 0.5, 1 2.5, 3.5 2.5)")
Lin1po1 <- gIntersection(Lin1, Po1)
coordinates(Lin1po1) 


plot(sl)
plot(l1, add=T, col="red", lwd=2)
plot(l2, add=T, col="blue", lwd=2)

gIntersects(l1, l2)
gTouches(l1, l2)

plot(createLine(x=cc[i:(i+1), 1], y=cc[i:(i+1), 2]))
