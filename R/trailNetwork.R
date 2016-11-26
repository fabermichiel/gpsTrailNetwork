# 
# Create a network of tracks from different input tracks and lines
# Similar overlapping tracks on a real world path are averaged and combined.
# Intersections are identified. Paths are created between those intersections, resulting in a network of trails.
# This network can be used to create maps, or for routing.

# Set up environment
rm(list=ls())

setwd("~/Documenten/scripts/trailNetwork_tmp")
source("functions.R")

loadPackage("rgdal")
loadPackage("rgeos")
loadPackage("sp")
loadPackage("tcltk")



# RD_New <- CRS("+init=epsg:28992")
# WGS84 <- CRS("+init=epsg:4326")
# WGS84.proj <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
#CRS(WGS84.proj)

# An import module of some kind should be added

### create the list of edges
# edges <- list()
# edges <- addEdge(createLine(x = c(1, 4, 6, 8, 8.5, 9, 13, 15, 16, 17, 19), 
#                             y = c(2, 2, 4, 1, 2.0, 5, 8 , 6 , 4 , 2 , 1), sort=TRUE), edges)
# edges <- addEdge(createLine(x = c(1, 2, 6, 7, 7.5, 8, 10 , 11, 12, 12.5, 14, 16, 18), 
#                             y = c(5, 4, 2, 3, 9.0, 7, 5 , 4 , 7 , 8 , 9 , 2 , 3), sort=TRUE), edges)
# edges <- addEdge(createLine(x = c(10, 1), 
#                             y = c(10, 1), sort=TRUE), edges)


# load gpx from disk
edges <- loadFiles()

####
# edges <- addEdge(createLine(x = c(5, 10), 
#                             y = c(9, 10), sort=TRUE))
# edges <- addEdge(createLine(x = c(9, 11), 
#                             y = c(10, 10), sort=TRUE))
# 


#### 
# track1 <- readOGR(dsn = "tracks/track1.gpx", layer = "route_points")
# edges <- addEdge(l = createLine(x = getCoordinatesOfPoints(track1), epsg = 4326), e = edges)
# 
# track2 <- readOGR(dsn = "tracks/track2.gpx", layer = "route_points")
# edges <- addEdge(l = createLine(x = getCoordinatesOfPoints(track2), epsg = 4326), e = edges)
# 
# track3 <- readOGR(dsn = "tracks/track3.gpx", layer = "route_points")
# edges <- addEdge(l = createLine(x = getCoordinatesOfPoints(track3), epsg = 4326), e = edges)
# 
# 


# transform to projected coordinates (default is 28992 = dutch national grid)
edges <- lapply(edges, FUN=transformToProjected)

beginEdges <- edges
#edges <- beginEdges


# plot(edges[[1]], type="b", col="green")#, xlim=c(6.11,6.14), ylim=c(52.42,52.45))#, xlim=c(40,60))
# plot(edges[[2]], add=T, type="b", col="green")
# plot(edges[[3]], add=T, type="b", col="green")
# # plot(edges[[4]], add=T, type="b", col="yellow")
# plot(edges[[5]], add=T, type="b", col="purple")
# 
# edges <- removeEdge(4)
# edges <- removeEdge(5)


t1 <- edges[[1]]
t2 <- edges[[2]]



#plot(t1, type="n", col="green")#, xlim=c(204000,206000), ylim=c(493000,496000))#, xlim=c(40,60)) # plots nothing, setup of plot
# plot(t1, type="b", col="green", add=T)
# plot(t2, type="b", col="green", add=T)
# plot(t3, type="b", col="green", add=T)

l1 <- createLine(head(getCoordinatesOfALine(t1)[c(1, 3, 4),]))
l2 <- createLine(head(getCoordinatesOfALine(t2)[c(1,2, 4), ]))

plot(l1, col="blue")
plot(l2, add=T, col="red")

edges[[1]] <- l1
edges[[2]] <- l2

#plotNumberOfPointsInRadius(e = edges, r = 2)

#Sys.sleep(0.5)

# We need to find the intersection between each edge. The resulting set of edges only intersects on incident nodes. 
# If this condition is met, we can reduce (average) two edges or a face sharing the same end points into one new edge.

plotEdges(edges)

edges <- reduceIntersections(e = edges)
plotEdges(edges)

edges <- cleanEdgeList(e = edges)
plotEdges(edges)
edges <- parallelReduction(e = edges) # needs a tolerance of 10 meter
plotEdges(edges)

edges <- cleanEdgeList(e = edges)
plotEdges(edges)

#lapply(edges, plotEdge, color="black")


#lapply(edges, plotEdgesBuffer, b=1.5)


#plot(edges[[6]], type="n", col="black", xlim=c(0,15), ylim=c(0,15))
# 
# for (i in 1:length(edges)) {
#   edge <- edges[i]
#   if (is.na(edges[i])) {
#     next 
#   } else {
#       plot(edges[[i]], add=T, type="b", col="red", lwd=4)
# #    text(getCoordinatesOfPoints(getStartPointOfLine(edges[[i]])), labels=i, pos=1)
#     #  Sys.sleep(0.1)
#   }
# }

edges <- serialReduction(e = edges)
edges <- cleanEdgeList(e = edges)

edges <- edgeContraction(tol = 2, e = edges)
#edges <- serialReduction(e = edges)
edges <- cleanEdgeList(e = edges)

#lapply(edges, plotEdge, color="red")


test <- function(e) {
  
  # 6 is begin van graph  -> naar begin 14
  # 9 ook, via 8 - > naar begin 14
  # van 7 naar 8/12 is inkomend. 12 -> begin 15
  
  indexOfLine <- c(6)
  indexOfOtherLine <- c(9,8)
  
  line <- e[[6]] 
  otherLine <- combineLines(e[[9]], e[[8]])
  
  #gDistance(line, otherLine, hausdorff = T) < 10
  if (gWithinDistance(line, otherLine, dist = 10, hausdorff = T)) {
    endPointOfLine <- getEndPointOfLine(line)
    endPointOtherLine <- getEndPointOfLine(otherLine)
    
    if (gWithin(endPointOfLine, endPointOtherLine) & gContains(endPointOfLine, endPointOfLine)) {
      # hopla, same endpoint
      
      newLine <- averageLines(l1 = line, l2 = otherLine)
      # find points on newline closest to start/endpoints of disconnected edges
      
      indexFirst <- getIndicesOfEqualPoints(p = getEndPointOfLine(l = e[[9]]), pos = "start", e = e)
      indexFirst <- indexFirst[!indexFirst %in% indexOfOtherLine]
      
      indexLast <- getIndicesOfEqualPoints(p = getEndPointOfLine(l = e[[9]]), pos = "end", e = e)
      indexLast <- indexLast[!indexLast %in% indexOfOtherLine]
      
      newPoint <- getClosestPointOnLine(p = getStartPointOfLine(e[[indexFirst]]), l = newLine)
      
      #cut newLine with these points and add points to disconnected edges
      newLines <- splitLineOnIntersection(l = newLine, p = getCoordinatesOfPoints(newPoint))
      
      e <- addEdge(l = newLines[[1]], e = e)
      e <- addEdge(l = newLines[[2]], e = e)
      
#       plot(newLines[[1]], add=T, type="b", col="red", lwd=4)
#       plot(newLines[[2]], add=T, type="b", col="purple", lwd=4)
#       plot(newPoint, add=T, pch=19, col="blue")
      
      

      e[[indexFirst]] <- changePointOfLine(l = e[[indexFirst]], p = newPoint, pos = "start")
      e[[indexLast]] <- changePointOfLine(l = e[[indexLast]], p = newPoint, pos = "end")
      
      
#       plot(e[[indexLast]], add=T, type="b", col="orange", lwd=2)
#       plot(e[[indexFirst]], add=T, type="b", col="orange", lwd=2)
      
      # remove all averaged edges
      for (i in indexOfLine) {
        e <- removeEdge(i = i, e=e)
      }
      for (i in indexOfOtherLine) {
        e <- removeEdge(i = i, e=e)
      }
      # even more
      
      # repeat
      # findGraph()
      # e <- test(e)
    } else {
      # continu searching with this line
    }
    
  } else {
    #stop searching because line is to far apart
    # continu with other line
  }
  
  # do next edge
  
  return(e)
}

newEdges <- test(e = edges)
newEdges <- serialReduction(e = newEdges)
newEdges <- parallelReduction(e = newEdges)
newEdges <- cleanEdgeList(e = newEdges)
lapply(newEdges, plotEdge, color="red")

# 
# for (i in 1:length(edges)) {
#   edge <- edges[i]
#   if (is.na(edges[i])) {
#     next 
#   } else {
# #    plot(edges[[i]], add=T, type="b", col="black", lwd=2)
#     plot(getStartPointOfLine(edges[[i]]), add=T, col=getVertexDegree(getStartPointOfLine(edges[[i]]), edges), pch=19)
#     text(getCoordinatesOfPoints(getStartPointOfLine(edges[[i]])), labels=i, pos=1)
#     #  Sys.sleep(0.1)
#   }
# }




# 
# plot(getFirstEdge(), add=T, lwd=3, col="red")
