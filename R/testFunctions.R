# testFunctions
rm(list=ls())

library(sp)
library(rgeos)
library(rgdal)
library(igraph)

setwd("~/Documenten/scripts/trailNetwork_tmp")
source("functions.R")

p <- createPoint(x=c(5.0, 5.5, 6.0), y=c(52, 52.5, 53))
class(p) == "SpatialPoints"
coordinates(p) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(p) == CRS("+init=epsg:28992")@projargs

getCoordinatesOfPoints(p) == as.data.frame(coordinates(p))

l <- createLine(x=c(5.0, 5.5, 6.0), y=c(52, 52.5, 53))
class(l) == "SpatialLines"
coordinates(l@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(l) == CRS("+init=epsg:28992")@projargs

getCoordinatesOfALine(l) == as.data.frame(coordinates(l))

sp <- getStartPointOfLine(l)
class(sp) == "SpatialPoints"
coordinates(sp) == coordinates(data.frame(x=c(5), y=c(52)))
proj4string(sp) == CRS("+init=epsg:28992")@projargs

ep <- getEndPointOfLine(l)
class(ep) == "SpatialPoints"
coordinates(ep) == coordinates(data.frame(x=c(6), y=c(53)))
proj4string(ep) == CRS("+init=epsg:28992")@projargs

rl <- reverseLine(l)
class(rl) == "SpatialLines"
coordinates(rl@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(6, 5.5, 5), y=c(53, 52.5, 52)))
proj4string(rl) == CRS("+init=epsg:28992")@projargs

sl <- sortLine(rl)
class(sl) == "SpatialLines"
coordinates(sl@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(sl) == CRS("+init=epsg:28992")@projargs

sl2 <- sortLine(l)
class(sl2) == "SpatialLines"
coordinates(sl2@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(sl2) == CRS("+init=epsg:28992")@projargs

cpl <- changePointOfLine(l, p=createPoint(4, 51), "s")
class(cpl) == "SpatialLines"
coordinates(cpl@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(4, 5.5, 6), y=c(51, 52.5, 53)))
proj4string(cpl) == CRS("+init=epsg:28992")@projargs

# combineer twee aan elkaar liggende lijnen l1 en l2. Beginpunt van l2 is eindpunt van l1
cl <- combineLines(l1 = l, l2 = rl)
class(cl) == "SpatialLines"
coordinates(cl@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6, 5.5, 5), y=c(52, 52.5, 53, 52.5, 52)))
proj4string(cl) == CRS("+init=epsg:28992")@projargs

#

# Split line l on intersection coordinates p
pi1 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=5, y=52), getCoordinatesOfALine(l))
pi1 == 1
spl1 <- splitLineOnIntersection(l, p=data.frame(x=5, y=52))
is.na(spl1[[1]])
class(spl1[[2]]) == "SpatialLines"
coordinates(spl1[[2]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(spl1[[2]]) == CRS("+init=epsg:28992")@projargs
png("test/spl1.png")
plot(l, type="b", col="green", lwd=5)
plot(createPoint(data.frame(x=5, y=52)), add=T, col="black", pch=19, cex=2)
#plot(spl1[[1]], add=T, col="blue", lwd=3)
plot(spl1[[2]], add=T, col="red", lwd=3)
dev.off()

pi2 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=5.25, y=52.25), getCoordinatesOfALine(l))
pi2 == 1
spl2 <- splitLineOnIntersection(l, p=data.frame(x=5.25, y=52.25))
class(spl2[[1]]) == "SpatialLines"
coordinates(spl2[[1]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.25), y=c(52, 52.25)))
proj4string(spl2[[1]]) == CRS("+init=epsg:28992")@projargs
class(spl2[[2]]) == "SpatialLines"
coordinates(spl2[[2]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5.25, 5.5, 6), y=c(52.25, 52.5, 53)))
proj4string(spl2[[2]]) == CRS("+init=epsg:28992")@projargs
png("test/spl2.png")
plot(l, type="b", col="green", lwd=5)
plot(createPoint(data.frame(x=5.25, y=52.25)), add=T, col="black", pch=19, cex=2)
plot(spl2[[1]], add=T, col="blue", lwd=3)
plot(spl2[[2]], add=T, col="red", lwd=3)
dev.off()

pi3 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=5.5, y=52.5), l = getCoordinatesOfALine(l))
pi3 == 2
spl3 <- splitLineOnIntersection(l, p=data.frame(x=5.5, y=52.5))
class(spl3[[1]]) == "SpatialLines"
coordinates(spl3[[1]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5), y=c(52, 52.5)))
proj4string(spl3[[1]]) == CRS("+init=epsg:28992")@projargs
class(spl2[[2]]) == "SpatialLines"
coordinates(spl3[[2]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5.5, 6), y=c(52.5, 53)))
proj4string(spl3[[2]]) == CRS("+init=epsg:28992")@projargs
png("test/spl3.png")
plot(l, type="b", col="green", lwd=5)
plot(createPoint(data.frame(x=5.5, y=52.5)), add=T, col="black", pch=19, cex=2)
plot(spl3[[1]], add=T, col="blue", lwd=3)
plot(spl3[[2]], add=T, col="red", lwd=3)
dev.off()

pi4 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=5.75, y=52.75), getCoordinatesOfALine(l))
pi4 ==2 
spl4 <- splitLineOnIntersection(l, p=data.frame(x=5.75, y=52.75))
class(spl4[[1]]) == "SpatialLines"
coordinates(spl4[[1]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 5.75), y=c(52, 52.5, 52.75)))
proj4string(spl4[[1]]) == CRS("+init=epsg:28992")@projargs
class(spl4[[2]]) == "SpatialLines"
coordinates(spl4[[2]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5.75, 6), y=c(52.75, 53)))
proj4string(spl4[[2]]) == CRS("+init=epsg:28992")@projargs
png("test/spl4.png")
plot(l, type="b", col="green", lwd=5)
plot(createPoint(data.frame(5.75, 52.75)), add=T, col="black", pch=19, cex=2)
plot(spl4[[1]], add=T, col="blue", lwd=3)
plot(spl4[[2]], add=T, col="red", lwd=3)
dev.off()

pi5 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=6, y=53), l = getCoordinatesOfALine(l))
pi5 == 3
spl5 <- splitLineOnIntersection(l, p=data.frame(x=6, y=53))
class(spl5[[1]]) == "SpatialLines"
coordinates(spl5[[1]]@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6), y=c(52, 52.5, 53)))
proj4string(spl5[[1]]) == CRS("+init=epsg:28992")@projargs
is.na(spl5[[2]])
png("test/spl5.png")
plot(l, type="b", col="green", lwd=5)
plot(createPoint(data.frame(x=6, y=53)), add=T, col="black", pch=19, cex=2)
plot(spl5[[1]], add=T, col="blue", lwd=3)
#plot(spl5[[2]], add=T, col="red", lwd=3)
dev.off()

# what happens if line goed back in x and y direction?
l1 <- createLine(x=c(5.0, 5.5, 6.0, 6.5, 7.0, 6.8, 6.8, 6.5, 6.0, 5.3), 
                 y=c(52.0, 52.5, 53.0, 53.1, 52.8, 52.5, 52.3, 52.3, 52.0, 51.8))

pi6 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=6.25, y=52.15), l = getCoordinatesOfALine(l1))
pi6 == 8
spl6 <- splitLineOnIntersection(l1, p=data.frame(x=6.25, y=52.15))
rl1 <- spl6[[1]]
rl2 <- spl6[[2]]
class(rl1) == "SpatialLines"
coordinates(rl1@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5, 5.5, 6, 6.5, 7.0, 6.8, 6.8, 6.5, 6.25), y=c(52, 52.5, 53, 53.1, 52.8, 52.5, 52.3, 52.3, 52.15)))
proj4string(rl1) == CRS("+init=epsg:28992")@projargs
class(rl2) == "SpatialLines"
coordinates(rl2@lines[[1]]@Lines[[1]]) == coordinates(data.frame(x=c(5.3, 6, 6.25), y=c(51.8, 52.0, 52.15)))
proj4string(rl2) == CRS("+init=epsg:28992")@projargs

png("test/spl6.png")
plot(l1, type="b", col="green", lwd=5)
plot(createPoint(data.frame(x=6.25, y=52.15)), add=T, col="black", pch=19, cex=2)
plot(spl6[[1]], add=T, col="blue", lwd=3)
plot(spl6[[2]], add=T, col="red", lwd=3)
dev.off()

# 
# gives expected error
# pi7 <- getPointIdOfLastPointOfIntersection(p=data.frame(x=5.28, y=51.7), l = getCoordinatesOfALine(l1))

#################
temp <- function(p, l) {
  d <- gDistance(p, l)
  
  plot(gBuffer(l, width = d+d/10), border="white")
  plot(l, type="b", add=T)
  plot(p, add=T)
  
  buff <- gBuffer(spgeom = p, width = d, quadsegs = 100)
  plot(buff, add=T)
  
  np <- getClosestPointOnLine(p, l)
  plot(np, add=T)
  nl <- createLine(rbind(getCoordinatesOfPoints(p), getCoordinatesOfPoints(np)))
  
  plot(nl, col="red",add=T)
  
  return(np)
}



########################################33
## average two lines
# simple test
l1 <- createLine(x=c(5.0, 5.5, 6.0, 6.5, 7.0), y=c(52, 52.5, 53.0, 53.1, 52.8))
l2 <- createLine(x=c(5.0, 5.2, 5.2, 5.3, 5.8, 6.2, 6.5, 7.0), y=c(52, 52.3, 52.6, 52.9, 53.2, 53.5, 53.6, 52.8))
la <- averageLines(l1, l2)
png("test/averageLine1.png")
plot(l2, type="n")
plot(l1, add=T, type="b", col="red", lwd=2)
plot(l2, add=T, type="b", col="blue", lwd=2)
plot(la, add=T, type="b", col="green", lwd=2)
dev.off()


# Outside with more points
l1 <- createLine(x=c(5.0, 5.5, 6.0, 6.5, 7.0, 6.8, 6.8, 6.5, 6.0, 5.3), 
                 y=c(52.0, 52.5, 53.0, 53.1, 52.8, 52.5, 52.3, 52.3, 52.0, 51.8))
l2 <- createLine(x=c(5.00, 5.20, 5.20, 5.30, 5.80, 6.20, 6.50, 7.20, 7.20, 7.00, 6.80, 6.50, 6.40, 6.20, 5.70, 5.50, 5.30), 
                 y=c(52.0, 52.3, 52.6, 52.9, 53.2, 53.5, 53.6, 52.8, 52.6, 52.5, 52.2, 52.1, 52.1, 51.7, 51.6, 51.7, 51.8))
la <- averageLines(l1, l2)
png("test/averageLine2.png")
plot(l2, type="n")
plot(l1, add=T, type="b", col="red", lwd=2)
plot(l2, add=T, type="b", col="blue", lwd=2)
plot(la, add=T, type="b", col="green", lwd=2)
dev.off()

# inside with more points.
# If closest point of outer line is further away then closest point at 'other side', the other side point is taken. This is incorrect.
# both path should be closer together then treshold dist.
l1 <- createLine(x=c(5.00, 5.12, 5.25, 5.37, 5.5, 5.67, 5.75, 5.87, 6.00, 6.12, 6.25, 6.37, 6.5, 6.67, 6.75, 6.87, 7.00, 6.90, 6.75, 6.70, 6.60, 6.50, 6.37, 6.25, 6.12, 6.00, 5.82, 5.65, 5.47, 5.30), 
                 y=c(52.0, 52.12, 52.25, 52.37, 52.4, 52.5, 52.75, 52.87, 53.0, 53.02, 53.05, 53.07, 53.1, 53.0, 52.9, 52.8, 52.75, 52.7, 52.6, 52.5, 52.5, 52.4, 52.3, 52.3, 52.2, 52.1, 52.0, 51.9, 51.9, 51.8))
l2 <- createLine(x=c(5.00, 5.20, 5.20, 5.30, 5.80, 6.20, 6.50, 7.20, 7.20, 7.00, 6.80, 6.50, 6.40, 6.20, 5.70, 5.50, 5.30), 
                 y=c(52.0, 52.3, 52.6, 52.9, 53.2, 53.5, 53.6, 52.8, 52.6, 52.5, 52.2, 52.1, 52.1, 51.7, 51.6, 51.7, 51.8))
la <- averageLines(l1, l2)
png("test/averageLine3.png")
plot(l2, type="n")
plot(l1, add=T, type="b", col="red", lwd=2)
plot(l2, add=T, type="b", col="blue", lwd=2)
plot(la, add=T, type="b", col="green", lwd=2)
dev.off()

e <- createEdgeList()
e <- addEdge(createLine(x=c(5.0, 5.5, 6.0), y=c(52, 52.5, 53)), e)
e <- addEdge(createLine(x=c(6.0, 6.1, 6.2), y=c(53, 52.5, 52)), e)

plot(combineLines(l1 = e[[1]], l2 = e[[2]]), type="n")
plot(e[[1]], add=T, type="b", col="red", lwd=5)
plot(e[[2]], add=T, type="b", col="green", lwd=5)

e <- serialReduction(e)
plot(getFirstEdge(e = e), add=T, type="b", col="black", lwd=2)
