
loadPackage <- function(package) {
  result <- require(package, character.only = T)
  if (!result) {
    install.packages(package)
    result <- require(package, character.only = T)
    if (!result) {
      error("Can't load and install package")
    }
  }
}

showEpsgs <- function(x=NULL) {
  names <- c("WGS 84", "RD New") 
  epsg <- c(4326, 28992)
  if (is.null(x)) {
    return(data.frame(names, epsg))
  }
}

getEpsg <- function(sp) {
  pj4 <- proj4string(sp)
  epsg <- as.integer(gsub(pattern = "epsg:", "", x=unlist(strsplit(unlist(strsplit(x = pj4, split = " +"))[1], "="))[2]))
  return(epsg)
}

createPoint <- function(x, y=NULL, epsg = 28992) {
  # input is vector with x coordinates and a vector with y coordinates OR x as dataframe with first column x coord and second column y coord
  # is x is a data.frame, a two column dataframe is assumed with first column of x coordinates en second column with y coordinates
  # epsg is EPSG of projection used. Default is dutch national grid
  # use a projected grid (units in meters or miles) and not a sperical (degrees) because of distance calculations
  # output is a spatialPoint
  # if epsg == NULL or NA, no projection will be used. This is the default of SpatialPoints()
  if (is.null(epsg)) {
    pj4 <- CRS(as.character(NA))
  } else {
    if (is.na(epsg)) {
      pj4 <- CRS(as.character(NA))
    } else {
      pj4 = CRS(paste("+init=epsg", epsg, sep=":"))
    }
  }
  if (is.data.frame(x)) {
    newPoint <- SpatialPoints(x, proj4string = pj4)
  } else {
    newPoint <- SpatialPoints(coords = data.frame(x=x, y=y), proj4string = pj4)
  }
  return(newPoint)
}

createLine <- function(x, y, doSort=TRUE, epsg = 28992) {
  # input is vector with x coordinates and a vector with y coordinates
  # is x is a data.frame, a two column dataframe is assumed with first column of x coordinates en second column with y coordinates
  # epsg is EPSG of projection used. Default is dutch national grid
  # use a projected grid (units in meters or miles) and not a sperical (degrees) because of distance calculations
  # output is a spatialLine
  pj4 = CRS(paste("+init=epsg", epsg, sep=":"))
  if (is.data.frame(x)) {
    newLine <- SpatialLines(list(Lines(list(Line(x)),ID=1)), proj4string = pj4)
  } else {
    newLine <- SpatialLines(list(Lines(list(Line(data.frame(x, y, stringsAsFactors = F))), ID=1)), proj4string = pj4)    
  }
  if (doSort) {
    newLine <- sortLine(newLine)
  }
  return(newLine)
}

transformToProjected <- function(sp, epsg = 28992) {
  res <- spTransform(sp, CRSobj = CRS(paste("+init=epsg", epsg, sep=":")))
  return(res)
}

getCoordinatesOfPoints <- function(p, asDataFrame=TRUE) {
  # input is a SpatialPoint
  # output is a dataframe when asDataFrame is True, otherwise a matrix with coordinatescoordinates <-p@coords
  coords <- coordinates(p)
  if (asDataFrame) {
    coords <- as.data.frame(coords, stringsAsFactors = F)
  }
  row.names(coords) <- 1:nrow(coords)
  return(coords)
}

getCoordinatesOfALine <- function(l, asDataFrame=TRUE) {
  # input is a SpatialLine
  # output is a dataframe when asDataFrame is True, otherwise a matrix with coordinates
  coords <-coordinates(l)
  if (asDataFrame) {
    coords <- as.data.frame(coords, stringsAsFactors = F)
  }
  colnames(coords) <- c("x","y")
  return(coords)
}

getCoordinates <- function(sp) {
  result <- NULL
  if (class(sp) == "SpatialPoints") {
    result <- getCoordinatesOfPoints(sp)
  } else {
    if (class(sp) == "SpatialLines") {
      result <- getCoordinatesOfALine(sp)
    } else {
      warning("Not a supported spatial object as input")
    }
  }
  return(result)
}

getStartPointOfLine <- function(l) {
  coordinatesOfLine <- getCoordinatesOfALine(l)
  newPoint <- createPoint(x = coordinatesOfLine$x[1], 
                          y = coordinatesOfLine$y[1])
  return(newPoint)
}

getEndPointOfLine <- function(l) {
  coordinatesOfLine <- getCoordinatesOfALine(l)
  numberOfPoints <- nrow(coordinatesOfLine)
  newPoint <- createPoint(x = coordinatesOfLine$x[numberOfPoints], 
                          y = coordinatesOfLine$y[numberOfPoints])
  return(newPoint)
}

reverseLine <- function(l) {
  coordsOfLine <- getCoordinatesOfALine(l)
  numberOfPoints <- nrow(coordsOfLine)  
  newLine <- createLine(x = coordsOfLine$x[numberOfPoints:1], 
                        y = coordsOfLine$y[numberOfPoints:1], 
                        doSort = FALSE, 
                        epsg = getEpsg(l))
  return(newLine)
}

sortLine <- function(l) {
  coordsOfLine <- getCoordinatesOfALine(l)
  numberOfPoints <- nrow(coordsOfLine)
  
  xOfFirstPoint <- coordsOfLine$x[1]
  xOfLastPoint <- coordsOfLine$x[numberOfPoints]
  if (xOfFirstPoint > xOfLastPoint) {
    newLine <- reverseLine(l)
  } else {
    newLine <- l
  }
  return(newLine)
}

# interaction with SpatialObjects 

changePointOfLine <- function(l, p, pos) {
  # changes first or last point of line with supplied point
  # input is SpatialLine, SpatialPoint and indication of sart or endpoint should be replaced
  coordinatesOfLine <- getCoordinatesOfALine(l)
  coordinatesOfPoint <- getCoordinatesOfPoints(p)
  if (pos %in% c("s", "start", "b", "begin")) {
    coordinatesOfLine[1,] <- coordinatesOfPoint 
  }
  if (pos %in% c("l", "last", "e", "end")) {
    coordinatesOfLine[nrow(coordinatesOfLine),] <- coordinatesOfPoint 
  }
  newLine <- createLine(coordinatesOfLine, doSort = FALSE)
  return(newLine)
}

combineLines <- function(l1, l2) {
  # input: (SpatialLines, SpatialLines)
  # output: SpatialLines
  coordinatesOfLine1 <- getCoordinatesOfALine(l1)
  coordinatesOfLine2 <- getCoordinatesOfALine(l2)
  
  coordinatesOfNewLine <- rbind(coordinatesOfLine1, coordinatesOfLine2[2:nrow(coordinatesOfLine2),])
  newLine <- createLine(coordinatesOfNewLine)
  return(newLine)
}

getPointIdOfLastPointOfIntersection <- function(p, l) {
  # input is a point and a line. Point is an intersection on line
  # find the rowid of the last point of a line (l) where point (p) intersects
  # ***both p and l are given as a dataframe with only x and y***
  # output is the row id of the last point before the intersection
  
  #   isEqual <- p$x == l$x & p$y == l$y
  #   if (any((isEqual))) {
  #     # punt ligt op node
  #     result <- which(isEqual)
  #   } else {
  #     lp1 <- which(p$x > l$x)
  #     lp3 <- which(p$y > l$y)
  #     
  #     lp2 <- which(p$x < l$x)
  #     lp4 <- which(p$y < l$y)  
  #     
  #     lp13 <- unique(lp1, lp3)
  #     lp24 <- unique(lp2, lp4)
  #     if (length(lp13) == 0 || length(lp24) == 0) {
  #       result <- NULL
  #     } else {
  #       result <- min(lp13[(lp13+1) %in% lp24])
  #     }
  #   }
  
  #  startPoint <- l[1,]
  
  #  d <- gDistance(createPoint(startPoint), createPoint(p))
  
  if (class(p)=="SpatialPoints") {
    p <- getCoordinatesOfPoints(p)
  }
    
  if (class(l)=="SpatialLines") {
    l <- getCoordinatesOfALine(l)
  }
  
  result <- NULL
  isEqual <- p$x == l$x & p$y == l$y
  if (any(isEqual)) {
    # punt ligt op node
    result <- which(isEqual)
  } else {
    
    numberOfPointsOnLine <- nrow(l)
    for (i in 1:(numberOfPointsOnLine-1)) {
      
      # loop gaat door tot eind van rij
      # hier wordt getest met i+1 wat een out of array gaat geven
      # echter code komt hier in principe nooit. Vang toch af 
      
      if (i == numberOfPointsOnLine) {
        # should never come to this point
        stop("point is not in line")       
      } else {
               
        # x en y moeten tussen x en y van i en i+1 liggen
        # Vind daarom eerst de hoogste (en dus ook laagste) x en y en kijk of p daar tussen ligt
        if (l$x[i] <= l$x[i+1]) {
          x1 <- l$x[i]
          x2 <- l$x[i+1]
        } else {
          x2 <- l$x[i]
          x1 <- l$x[i+1]
        }
        # Is this what i want?
        # I switch y, but i could be that x is not switched...
        if (l$y[i] <= l$y[i+1]) {
          y1 <- l$y[i]
          y2 <- l$y[i+1]
        } else {
          y2 <- l$y[i]
          y1 <- l$y[i+1]
        }
        if (x1 <= p$x & x2 >= p$x & y1 <= p$y & y2 >= p$y) {
          # punt ligt tussen nodes
          
          # is het ook het goede punt?
          # test dit door te kijken of gezocht punt ook op de lijn ligt
          newLine <- createLine(rbind(l[i,],l[i+1,]))
          #          lines(newLine, col="yellow", lwd=3)
          doIntersect <-  gIntersects(newLine, gBuffer(createPoint(p), width = 0.0003))
          if (doIntersect) {
            # correct point
            result <- i
            break
          } else {
            next
          }
        }
        #      }
      } 
    }
  }
  return(result)
}

splitLineOnIntersection <- function(l, p, doSort=TRUE) {
  
  # input is the SpatialLine (l) to split and the point (p) where to split the line
  # Point is SpatialPoint or dataframe with coordinates
  # output is a list of two SpatialLines split on the point 
  
  # Example
  # l <- createLine(x=c(5.0, 5.5, 6.0), y=c(52, 52.5, 53))
  # p <- createPoint(x=5.1, y=52.1)
  # res <- splitLineOnIntersection(l, p) 
  #   
  # plot(l, type="b")
  # plot(p, pch=19, col="red", ,cex=2, add=T)
  # plot(res[[1]], col="red", lwd=2, add=T) 
  # plot(res[[2]], col="blue", lwd=2, add=T)
  
  
  coordsOfLine <- getCoordinates(l)
  if (class(p)=="SpatialPoints") {
    coordsOfPoint <- getCoordinates(p)
  } else {
    coordsOfPoint <- p
  }
  numberOfPointsOnLine <- nrow(coordsOfLine)
  
  idOfPoint <- getPointIdOfLastPointOfIntersection(p, l)
  if (is.null(idOfPoint)) {
    stop("No point found before intersecting point. Please check and debug")
  }
  if (idOfPoint == 1 & all(coordsOfLine[idOfPoint,] == coordsOfPoint)) {# intersect on startPoint
    nl1 <- NA
    nl2 <- l
  } else {
    if (idOfPoint != 1 & idOfPoint != numberOfPointsOnLine & all(coordsOfLine[idOfPoint, ] == coordsOfPoint)) {  # intersect on a middle node
      # split line equally
      nl1 <- createLine(x=c(coordsOfLine$x[1:idOfPoint]), 
                        y=c(coordsOfLine$y[1:idOfPoint]), doSort=doSort)  
      nl2 <- createLine(x=c(coordsOfLine$x[(idOfPoint):numberOfPointsOnLine]), 
                        y=c(coordsOfLine$y[(idOfPoint):numberOfPointsOnLine]), doSort=doSort)
    } else {
      if (idOfPoint == numberOfPointsOnLine & all(coordsOfLine[idOfPoint, ] == coordsOfPoint)) { # intersect on end point
        nl1 <- l
        nl2 <- NA
      } else {     
        # not on a node of the line
        nl1 <- createLine(x=c(coordsOfLine$x[1:idOfPoint], p$x), 
                          y=c(coordsOfLine$y[1:idOfPoint], p$y), doSort=doSort)  
        nl2 <- createLine(x=c(p$x, coordsOfLine$x[(idOfPoint+1):numberOfPointsOnLine]), 
                          y=c(p$y, coordsOfLine$y[(idOfPoint+1):numberOfPointsOnLine]), doSort=doSort)
      }
    }
  }
  return(list(nl1, nl2))  
}

cutLineOnIntersections <- function(sl, int, doSort=TRUE) {
  # Spatialline, data.frame -> list with SpatialLines
  
  # split a SpatialLine (sl) on the given coordinates (int)
  
  # example 
  
  # x1 <- c(1, 3,   5, 6, 7)
  # y1 <- c(1, 2.5, 2, 1, 1.5)
  # x2 <- c(8, 6, 5, 1, 0)
  # y2 <- c(3, 2, 1, 2, 2)
  
  # Is there also a easier way of doing this?
  # s1 <- createLine(x1,y1)
  # s2 <- createLine(x2,y2)
  
  # int <- gIntersection(s1, s2)
  
  # plot(  s1, type="b", col="red",  lwd=3, xlim=c(0,10))
  # lines( s2, type="b", col="blue", lwd=3)
  # points(int, pch=19,   col="black") 
  
  #   result <- cutLineOnIntersections(s1, int)
  #   
  #   lines(result[[1]], col="green", lwd=3)
  #   lines(result[[2]], col="green", lwd=6)
  #   lines(result[[3]], col="green", lwd=9)
  #   
  
  coordinatesOfLine <- getCoordinates(sl)
  # coordinatesOfIntersections <- getCoordinates(int) # int is already a data frame
  coordinatesOfIntersections <- int
  numberOfIntersections <- nrow(coordinatesOfIntersections)
  
  result <- list()
  for (i in numberOfIntersections:1) {
    if (i==numberOfIntersections) {
      tmp <- splitLineOnIntersection(l = sl, p = int[i,], doSort = doSort)
      result[[i+1]] <- tmp[[2]]
    } else {
      tmp <- splitLineOnIntersection(tmp[[1]], int[i,], doSort = doSort)
      result[i+1] <- tmp[[2]]
    }
    if (i==1) {
      result[1] <- tmp[[1]]    
    }
  }  
  return(result)
}

splitIntersectingLines <- function(s1, s2) {
  # Spatialline, SpatialLine -> list with SpatialLines
  
  # split two (intersecting) SpatialLines (sl1 and sl2) on the shared coordinates
  
  # example 
  
  # x1 <- c(1, 3,   5, 6, 7)
  # y1 <- c(1, 2.5, 2, 1, 1.5)
  # x2 <- c(8, 6, 5, 1, 0)
  # y2 <- c(3, 2, 1, 2, 2)
  
  # Is there also a easier way of doing this?
  # s1 <- createLine(x1,y1)
  # s2 <- createLine(x2,y2)
  
  # int <- gIntersection(s1, s2)
  
  # plot(  s1, col="red",  lwd=3, xlim=c(0,10))
  # lines( s2, col="blue", lwd=3)
  # points(int, pch=19,   col="black") 
  
  #   result <- splitIntersectingLines(s1, s2)
  #   
  #   lines(result[[1]], col="green", lwd=3)
  #   lines(result[[2]], col="green", lwd=6)
  #   lines(result[[3]], col="green", lwd=9)
  #   lines(result[[4]], col="orange", lwd=3)
  #   lines(result[[5]], col="orange", lwd=6)
  #   lines(result[[6]], col="orange", lwd=9)
  
  int <- gIntersection(s1, s2)
  result1 <- cutLineOnIntersections(s1, int)
  result2 <- cutLineOnIntersections(s2, int)
  
  result <- c(result1, result2)
  
}

# To Do
cutLineOnSelfIntersection <- function(sl) {
  # SpatialLine -> list of SpatialLines
  
  # Find all self intersections of a spatialLine. If found, split the line on the self intersection(s)
  # Return a Spatial(Multi)Line with the split segments 
  
  ## example:
  #   sl <- readWKT("LINESTRING(1 1, 3 3, 1 3, 3 1)")
  #   sll <- cutLineOnSelfIntersection(sl)
  #   sll_reference <- readWKT("MULTILINESTRING((1 1, 2 2),(2 2, 3 3, 1 3, 2 2),(2 2, 3 1))")
  #   gEquals(sll, sll_reference)
  #   plot(sll, type="b")
  
  cc <- getCoordinatesOfALine(sl)
  
  int <- matrix(FALSE, nrow = nrow(cc), ncol = nrow(cc))
  for (i in 1:(nrow(cc)-2)) {
    l1 <- SpatialLines(LinesList = list(Lines(slinelist = list(Line(coords = data.frame(x=cc[i:(i+1), 1], y=cc[i:(i+1), 2]))), ID = 1)))
    
    for (j in (i+1):(nrow(cc)-1)) {
      l2 <- SpatialLines(LinesList = list(Lines(slinelist = list(Line(coords = data.frame(x=cc[j:(j+1), 1], y=cc[j:(j+1), 2]))), ID = 2)))
      int[i, j] <- gIntersects(l1, l2) & !gTouches(l1, l2)
    }
  }
  int[11] <- TRUE
  ii <- which(int)
  for (n in 1:length(ii)) {
    i <- ii[n] %% nrow(cc)
    j <- ceiling(ii[n] / nrow(cc))   
  }
}

getDouble <- function(coords) {
  doubleIndex <- duplicated(coords) 
  if (any(doubleIndex)) {
    # get both indices
  } 
}

findSelfIntersections <- function(sl) {
  # Spatialline -> list with list of spatialpoint and indices of crossing point coordinates
  
  # find any self intersection of a spatialline
  
  # example
  # x <- c(1, 3,   5, 6, 7,   9, 8, 6, 5, 1, 0)
  # y <- c(1, 2.5, 2, 1, 1.5, 2, 3, 2, 1, 2, 2)
  # coords <- data.frame(x, y)
  
  # sl <- SpatialLines(list(Lines(list(Line(coords)), ID=1)))
  # plot(sl, type="b", col="green")
  
  # int <- findSelfIntersections(sl)
  # plot(int[[1]][[1]], pch=19, add=T, col="red")
  # plot(int[[2]][[1]], pch=19, add=T, col="red")
  
  #============================
  
  # get coordinates of the line
  coords <- getCoordinatesOfALine(sl)
  
  # for every line segment, check if it crosses with another line segment
  result <- list()
  #result <- NULL
  
  # To do
  # Find double coordinates
  
  for (p in 1:(nrow(coords)-2)) {
    # get the first line between two consecutive points
    l1 <- SpatialLines(list(Lines(list(Line(coords[p:(p+1),])),ID=1)))
    
    for (q in (p+1):(nrow(coords)-1)) { 
      # q <- p+1
      # get the second line from next coordinates
      l2 <- SpatialLines(list(Lines(list(Line(coords[q:(q+1), ])),ID=2)))
      
      # do both lines intersect
      intersection <- getIntersection(l1, l2)
      # To do:
      # the result should be a list of lists with spatialPoint, p and q
      
      
      if (!is.null(intersection)) {
        int <- list(intersection, p, q)
        
        result[[length(result)+1]] <-  int
        
        
        #cop <- getCoordinatesOfPoints(intersection)
        # int <- rbind(int, cbind(p, q, cop))  
      }    
    }
  }
  return(result)
}

getIntersection <- function(sl1, sl2) {
  # Spatialline, SpatialLine -> SpatialPoint
  
  # returns the intersection points of two spatiallines.
  # if none or not as SpatialPoints, then return NULL
  
  
  result <- NULL
  if (gIntersects(sl1, sl2) & !gTouches(sl1, sl2)) {
    int <- gIntersection(sl1, sl2)
    if (class(int) == "SpatialPoints") {
      result <- int
      
      # to get exactly same result as from createPoint(), create points from coordinates of intersection
      result <- createPoint(getCoordinates(result))
    }
  }
  return(result)
}

convertInt <- function(int) {
  #  matrix -> matrix 
  
  # prepare int for intersections
  
  # Example
  # example
  # x <- c(1, 3,   5, 6, 7,   9, 8, 6, 5, 1, 0)
  # y <- c(1, 2.5, 2, 1, 1.5, 2, 3, 2, 1, 2, 2)
  # coords <- data.frame(x, y)
  
  # sl <- SpatialLines(list(Lines(list(Line(coords)), ID=1)))
  # plot(sl, type="b", col="green")
  
  # int <- findSelfIntersections(sl)
  # newInt <- convertInt(int)
  
  #--------------
  
  newInt <- rbind(int[, c(1, 3, 4)], int[, c(2, 3, 4)])
  newInt <- newInt[order(newInt[,"p"]),]
  
  return(newInt)
}

isLoop <- function(sl) {
  # SpatialLine -> Boolean
  
  # Given a spatialline, this function checks if it is a loop by testing if startpoint is equal to endpoint
  
  # Example
  # x <- c(5.5, 6, 7,   9, 8, 6, 5.5)
  # y <- c(1.5, 1, 1.5, 2, 3, 2, 1.5)
  # coords <- data.frame(x, y)
  
  # sl <- SpatialLines(list(Lines(list(Line(coords)),ID=1)))
  # plot(sl, type="b", col="green")
  
  # isLoop(sl)
  
  #--------------
  
  loop <- FALSE
  coords <- as.data.frame(coordinates(sl))
  
  # get startpoint
  startpoint <- getStartPointOfLine(sl)
  
  # get endpoint
  endpoint <- getEndPointOfLine(sl)
  
  if (gEquals(startpoint, endpoint)) {
    loop <- TRUE
  }
  
  return(loop)
}

splitLoop <- function(sl, verbose=FALSE) {
  # Spatialline -> list with Spatiallines
  
  # Split the loop into two part. The split is made on the point furthest away from startpoint
  
  # example
  # x <- c(5.5, 6, 7,   9, 8, 6, 5.5)
  # y <- c(1.5, 1, 1.5, 2, 3, 2, 1.5)
  # coords <- data.frame(x, y)
  
  # sl <- SpatialLines(list(Lines(list(Line(coords)),ID=1)))
  # plot(sl, type="b", col="green")
  # points( as.data.frame(coordinates(sl))[1, ], pch=19)
  # result <- splitLoop(sl)
  # lines(result[[1]], col="red")
  # lines(result[[2]], col="blue")
  # points( as.data.frame(coordinates(sl))[result[[3]], ], pch=19)
  
  #============================
  
  result <- list()
  coords <- getCoordinates(sl)
  # get startpoint
  p1 <- getStartPointOfLine(sl)
  distance <- 0
  pointId <- NA
  for (p in 2:(nrow(coords)-2)) {
    # point to check
    p2 <- coords[p, ]
    dp <- gDistance(SpatialPoints(p1), SpatialPoints(p2))
    if (dp > distance) {
      # new distance
      distance <- dp
      # point furthest away from startpoint
      pointId <- p 
    }
  }
  
  # split line on point  
  l1 <- createLine(x = coords[1:pointId, ])
  l2 <- createLine(x = coords[pointId:nrow(coords), ])
  
  result[1] <- l1
  result[2] <- l2
  if (verbose) {
    result[3] <- pointId 
    result[4] <- distance
  }
  return(result)
}

getLineWithMostPoints <- function(l1, l2) {
  # input are two spatialLines
  # output is line_id with most points
  if (nrow(getCoordinatesOfALine(l1)) >= nrow(getCoordinatesOfALine(l2))) {
    return(1)
  } else {
    return(2)
  }
}

getClosestPointOnLine <- function(p, l) {
  #---------
  # if someone knows a better way to get closest point on line in R, please change
  #---------
  
  # Create a buffer from a point with a radius equal to shortest distance to line
  # intersect of buffer with line is target point
  # intersection comes up with two near equal points
  # average those two to get one point, error is small
  d <- gDistance(p, l)
  tempBuffer <- gBuffer(p, width=d+d/1000, quadsegs=1000)
  #plot(tempBuffer, add=T)
  
  newPoint <- gIntersection(tempBuffer, l)
  
  # added [,c("x","y")] because if a point is in the middle of a corner, there are two points closest to line. Only the first is used this way.
  newPoint <- createPoint(as.data.frame(t(colMeans(getCoordinatesOfPoints(p = newPoint)[,c("x","y")])), stringsAsFactors = F))
  #plot(newPoint, add=T) 
  return(newPoint)
}

getClosestTrackPointOnLine <- function(p, l) {
  coordinatesOfPoint <- getCoordinatesOfPoints(p)
  coordinatesOfLine <-  getCoordinatesOfALine(l)
  
  d <- sqrt((coordinatesOfLine$x-coordinatesOfPoint$x)^2+(coordinatesOfLine$y-coordinatesOfPoint$y)^2) # find all distances to current point 'x'
  selPointShort <- d == min(d) # search closest point
  #  d1 <- sum(coordinatesOfPoint$x,coordinatesOfLine$x[selPointShort])/2 # calculate average Y
  #  d2 <- sum(coordinatesOfPoint$y,coordinatesOfLine$y[selPointShort])/2 # calculate average X
  #  newPoint <- createPoint(d1, d2)
  newPoint <- createPoint(coordinatesOfLine[selPointShort,])
  
  #selP <- createPoint(coordinatesOfLine[selPointShort,])
  
  #plot(selP, add=T, pch=19, col="blue")
  #plot(np, add=T, pch=19, col="green")
  
  return(newPoint)
}

getPointsOnOtherLine <- function(l1, l2, fixToPoint=FALSE) {
  # input are two spatialLines: l1 is line with the most points
  # output is a dataframe with the set of pointcoords equal to size of points in line with most points. 
  # Output consist of pointcoords closets on l2 to points of l1
  # If fixToPoint is TRUE, only track points on l2 are used
  
  coordOfLine <- getCoordinatesOfALine(l1)
  numberOfPoints <- nrow(coordOfLine)
  
  # Set up a data frame for the points. First and last point are known, l1 and l2 have equal begin/end points
  newPoints <- as.data.frame(matrix(NA, nrow=numberOfPoints, ncol=2), stringsAsFactors = F)
  colnames(newPoints) <- c("x", "y")
  newPoints[1, ] <- coordOfLine[1,]
  newPoints[numberOfPoints, ] <- coordOfLine[numberOfPoints,]
  newPoint <- NULL
  
  for (pt in 2:(numberOfPoints-1)) {
    #pt <- pt+1
    print(pt)
    pointToLookFrom <- createPoint(coordOfLine[pt,])  
    # points(pointToLookFrom, col="green", pch=19)
    
    # When we search line 2 for closest points, we shouldn't look for points, behind the last found point. 
    # If we do, then points can be found after last point which result in a strange loop or bend, averaging twice with same points.
    # So l2 should be shortend/cut so only non-found linesegment/points are included in the search
    # But we should stop clipping if last found point is the endpoint
    if (is.null(newPoint)) {
      l2tmp <- l2
    } else { 
      if (!gEquals(getEndPointOfLine(l2), newPoint)) { # if not last point
        l2tmp <- cutLineOnIntersections(sl = l2tmp, int = getCoordinates(newPoint), doSort = FALSE)[[2]] # cut the line and use only the second/last line
      } else {
        # last point is only point left, thus last point is newPoint
        # fill all remaining slots with last found point and exit the for loop
        newPoints[pt:(numberOfPoints-1),] <- getCoordinates(newPoint)
        break
      }
    }
    
    # lines(l2tmp, col="orange", lwd=2*pt)
    
    if (fixToPoint) {
      # use only the track points from l2
      newPoint <- getClosestTrackPointOnLine(p = pointToLookFrom, l = l2tmp)
    } else {
      newPoint <-getClosestPointOnLine(p = pointToLookFrom, l = l2tmp)
    }    
    # points(newPoint, col="yellow", pch=19, cex=2*pt)
    
    newPoints[pt, ] <- getCoordinates(newPoint)

  }
  return(newPoints)
} 

averageLineWithPoints <- function(l, p) {
  # input is SpatialLine with most points and a dataframe with the set of pointcoords closest to line
  # output is new average SpatialLine
  coordsOfLine <- getCoordinatesOfALine(l)
  newLine <- createLine(x=rowMeans(cbind(p$x, coordsOfLine$x)), 
                        y=rowMeans(cbind(p$y, coordsOfLine$y)))
  return(newLine)
}

averageLines <- function(l1, l2, fixToPoint=FALSE) {
  # input are two spatialLines
  # output is spatialLine
  if (getLineWithMostPoints(l1, l2) == 1 ) {
    more <- l1
    less <- l2
  } else {
    more <- l2
    less <- l1
  }
  newPoints <- getPointsOnOtherLine(l1 = more, l2 = less, fixToPoint = fixToPoint)
  newLine <- averageLineWithPoints(more, newPoints)
  return(newLine)
}

# interaction with edges

createEdgeList <- function() {
  edges <- list()
  return(edges)
}

addEdge <- function(l, e) {
  # add a edge to the list of edges
  e[length(e)+1] <- l
  return(e)
}

removeEdge <- function(i, e) {
  # remove a edge
  e[i] <- NA
  return(e)
}

getFilledEdges <- function(e) {
  res <- which(!unlist(lapply(e, is.na)))
  return(res)
}

getEmptyEdges <- function(e) {
  res <- which(unlist(lapply(e, is.na)))
  return(res)
}

cleanEdge <- function(i, e) {
  # delete a edge from the list of edges
  e[i] <- NULL
  return(e)
}

cleanEdgeList <- function(e) {
  # removes all deleted edges from the list of edges
  e <- cleanEdge(getEmptyEdges(e), e)
  
  #   numberOfEdges <- length(e)
  #   for (edge in 1:numberOfEdges) {
  #     if (is.na(e[edge])) {
  #       e[edge] <- NULL
  #       e <- cleanEdgeList(e = e)
  #       return(e)
  #     }
  #   } 
  return(e)
}

getEdge <- function(i, e, na.allowed = FALSE) {
  index <- i
  if (!na.allowed & is.na(e[[index]])) {
    index <- index+1
    edge <- getEdge(index+1, e)
  } else {
    edge <- e[[index]]
  }
  return(edge)
}

getFirstEdge <- function(wantIndex = FALSE, e) {
  # getFirstEdge <- function(wantIndex = FALSE, e) {
  #   numberOfEdges <- length(e)
  #   for (edge in 1:numberOfEdges) {
  #     if (is.na(e[edge])) {
  #       next
  #     } else {
  #       newLine <- e[[edge]]
  #       index <- edge
  #       break
  #     }
  #   }
  #   if (wantIndex) {
  #     result <- index
  #   } else {
  #     result <- newLine
  #   }
  #   return(result)
  # }
  
  index <- 1
  newLine <- getEdge(index, e)
  if (is.na(newLine)) {
    index <- index+1
    newLine <- getEdge(index, e)
  }
  if (wantIndex) {
    result <- index
  } else {
    result <- newLine
  }
  return(result)
}

splitLine <- function(i, ps, e) {
  # i is index of line to split
  # ps is coordinate set of intersections
  # e is set of edges
  
  # split the line in n+ps parts
  for (n in 1:nrow(ps)) { # for the first row, get input line
    if (n==1) { # for the first intersection, split the input line in two and remove input from edges
      newLines <- splitLineOnIntersection(l=e[[i]], p=ps[n,])
      e <- removeEdge(i, e) # is netjes maar langzaam. Sneller is e[i] <- NA
    } else { # for all other lines, get the result of previous split
      newLines <- splitLineOnIntersection(l=newLines[[2]], p=ps[n,])
    }
    
    # add the new first split to the edges
    e <- addEdge(newLines[[1]], e) # is netjes maar langzaam. Sneller is e[length(e)+1] <- newLines[[1]]    
    
    if (n == nrow(ps)) { # if at the end, add remaining line segment
      e <- addEdge(newLines[[2]], e) # is netjes maar langzaam. Sneller is e[length(e)+1] <- newLines[[2]]
    }       
  }
  return(e)
}

getIndicesOfEqualPoints <- function(p, pos, e) {
  # input is a SpatialPoint to look for, an agrument for kind of point (start of end) and a set of edges
  # output is vector of index of edges with same point on give location
  numberOfEdges <- length(e)
  index <- NULL
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      if (pos %in% c("s", "start", "b", "begin")) {
        point <- getStartPointOfLine(e[[edge]])
      } else {
        if (pos %in% c("l", "last", "e", "end")) {
          point <- getEndPointOfLine(e[[edge]])  
        } else {
          stop("Give a correct argument what kind of point to look for: start or end")
        }
      }
      if (gEquals(p, point)) {
        index <- c(index, edge)
      }
    }
  }
  return(index)
}

getVertexDegree <- function(p, e) {
  numberOfEdges <- length(e)
  count <- 0
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      startPoint <- getStartPointOfLine(e[[edge]])
      if (gEquals(p, startPoint)) {
        count <- count + 1
      }
      endPoint <- getEndPointOfLine(e[[edge]])
      if (gEquals(p, endPoint)) {
        count <- count + 1
      }
    }
  }
  return(count)
}

getIndexOfEdgesWithVertexDegree <- function(d, e) {
  # returns the index numbers of edges with certain incident degree
  # only startpoints are evaluated
  # not found endpoints have degree of 1
  result <- unlist(lapply(lapply(e, getStartPointOfLine), getVertexDegree, e = e)) == d
  if (any(result)) {
    result <- which(result)
  } else {
    result <- NULL
  }
  return(result)
}

getNumberOfPointsInRadiusFromPoint <- function(p, e, r) {
  buffer <- gBuffer(p, width = r )
  numberOfEdges <- length(e)
  count <- 0
  
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      coordinatesOfLine <- getCoordinatesOfALine(e[[edge]])
      numberOfCoordinates <- nrow(coordinatesOfLine)
      for (point in 1:numberOfCoordinates) {
        if (gWithin(createPoint(coordinatesOfLine[point,]), buffer)) {
          count <- count + 1
        }
      }
    }
  }
  #  plot(buffer, add=T)
  
  return(count)
}

getLineWithStartPoint <- function(p, e) {
  # Input is spatialpoint and set of edges
  # output is index of SpatialLine with same startPoint as p
  # only the first line found is returned
  numberOfEdges <- length(e)
  result <- NULL
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      startPoint <- getStartPointOfLine(e[[edge]])
      if (gEquals(p, startPoint)) {
        result <- edge
      }
    }
  }
  return(result)
}

getLineWithStartPointLapply <- function(e, p) {
  # Does not work
  helperFunction <- function(l, p) {
    if (is.na(l)) {
      result <- NA
    } else {
      startPoint <- getStartPointOfLine(l)
      if (gEquals(p, startPoint)) {
        result <- 1
      }  else {
        result <- NA
      }
    }
    return(result)
  }
  res <- lapply(e, FUN=helperFunction, p=p)
  result <- min(which(!is.na(res)))
  return(result)
}

getIndexNextLine <- function(i, e) {
  indexOfNextLine <- getIndicesOfEqualPoints(p = getEndPointOfLine(l = e[[i]]), pos = "start", e = e)
  return(indexOfNextLine)
}

touchPointOnStartOrEndNode <- function(l1, l2, p) {
  
  coordinatesOfPoint <- getCoordinatesOfPoints(p)
  numberOfCoordsOfPoint <- nrow(coordinatesOfPoint)
  sp <- getCoordinatesOfPoints(getStartPointOfLine(l1))
  ep <- getCoordinatesOfPoints(getEndPointOfLine(l1))
  spo <- getCoordinatesOfPoints(getStartPointOfLine(l2))
  epo <- getCoordinatesOfPoints(getEndPointOfLine(l2))
  
  if (numberOfCoordsOfPoint == 1) {
    
    equalStartPoint <- all(sp == coordinatesOfPoint)
    equalEndPoint <- all(ep == coordinatesOfPoint)
    equalStartPointOther <- all(spo == coordinatesOfPoint)
    equalEndPointOther <- all(epo == coordinatesOfPoint)
    
    doTouch <- (equalStartPoint | equalEndPoint) & (equalStartPointOther | equalEndPointOther)
  } else {
    if (all(sp == coordinatesOfPoint[1,])) {
      doTouch <- FALSE 
    } else {
      doTouch <- TRUE
    }
  }
  return(doTouch)
}

findFace <- function(e) {
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      index <- NULL
      endPointsLine <- NULL
      line <- e[[edge]]
      startPoint <- getStartPointOfLine(l = e[[edge]])
      endPoint <- getEndPointOfLine(l = e[[edge]])
      dd <- data.frame(c(index=integer(), sp=list()))
      indexLine <- c(index, edge)
      endPointsLine <- c(endPointsLine, endPoint)
      
      indexOfOtherLine <- getIndicesOfEqualPoints(p = startPoint, "start", e = e)
      indexOfOtherLine <- indexOfOtherLine[!indexOfOtherLine %in% edge] # remove self from index
      lastPointOfOtherLine <- getEndPointOfLine(l = e[[indexOfOtherLine]])
      # This endpoint can't be the same as endpoint. If that is the case, then it would be reduced by paralell reduction
      
      # What should be done:
      # Get next line from line
      # add to line
      # check dist
      # check endpoint with endpoints of other line
      # check next line from second line
      # add to line 
      # check dist
      # compare endpoint with endpoints of line
      
      # iterate until found or dist is too large
      
      # Get next line from line
      # get indices for next lines where endPoint is startpoint
      indexOfNextLine <- getIndicesOfEqualPoints(p = endPoint, pos = "start", e = e)
      for (i in indexOfNextLine) {
        secondLine <- e[[i]]
        newLine <- combineLines(line, secondLine)
        endPointOfNewLine <- getEndPointOfLine(l = newLine)
        
        if (dist < tol) {
          # keep in mem
          # check ends
        }
        
      }
      
      
      
      indexOfOtherSecondLines <- getIndicesOfEqualPoints(p = lastPointOfOtherLine, "start", e)
      lastPointOfOtherSecondLine <- getEndPointOfLine(l = e[[indexOfOtherSecondLines]])
      
      indexOfOtherThirdLines <- getIndicesOfEqualPoints(p = lastPointOfOtherSecondLine, "start", e)
      
      for (n in 1:length(indexOfOtherThirdLines) ) {
        lastPointOfThirdOtherLine <- getEndPointOfLine(l = e[[n]])
        if (gEquals(endPoint, lastPointOfThirdOtherLine)) {
          # match
          # collect indices and create spatialLines
          indexOfOtherLines <- c(indexOfOtherLine, indexOfOtherSecondLines, indexOfOtherThirdLines[n])
          
          for (n in 1:(length(indexOfOtherLines)-1)) {
            if (n == 1) {
              newOtherLine <- combineLines(l1 = e[[indexOfOtherLines[n]]], e[[indexOfOtherLines[n+1]]])
            } else {
              newOtherLine <- combineLines(l1 = newOtherLine, e[[indexOfOtherLines[n+1]]])
            }
          }
          
          newLine <- e[[edge]]
          
          newAverageLine <- averageLines(l1 = newLine, l2 = newOtherLine)
          #plot(newAverageLine, add=T, col="blue")
          
          # delete the old lines
          # todo
          e <- removeEdge(edge, e)
          
        }
      }
      
    }
  }
  return()
}

createGraphEdges <- function(x) {
  # Convert the edgesInfo to a vector with nodes
  # edgesInfo is a dataframe with 'start' and 'end' point id's of the edges  
  # return a 'edgelist'
  e <- NULL
  for (i in 1:nrow(x))
    e <- c(e, x$start[i], x$end[i])
  return(e)
}

createGraph <- function(e) {
  # creates a graph object from an 'edgelist'
  g <- graph(e)
  return(g)
}

convertToGraph <- function(x) {
  # creates a graph object from edgesInfo
  # edgesInfo is a dataframe with 'start' and 'end' point id's of the edges  
  # returns a graphobject
  g <- graph.edgelist(as.matrix(edgesInfo[,c("startnode", "endnode")]))
  return(g)
}

getExtentOfEdges <- function(e) {
  b <- NULL
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      bp <- as.data.frame(t(e[[edge]]@bbox))
      
      maxX <- max(b[2, 1], bp[2, 1])
      maxY <- max(b[2, 2], bp[2, 2])
      minX <- min(b[1, 1], bp[1, 1])
      minY <- min(b[1, 2], bp[1, 2])
      
      b <- data.frame(x = c(minX, maxX), y =c(minY, maxY))    
    } 
  }
  return(b)
}


# plot

plotLines <- function(l1, l2, l3=NULL, l4=NULL) {
  # This function plots up to four SpatialLines with colors red, blue, green and purple
  
  # get the coorinates of the base lines
  cc1 <- getCoordinates(l1)
  if (!is.null(l2)) {
    cc2 <- getCoordinates(l2)  
  } else {
    cc2 <- data.frame(x=NA, y=NA)
  }
  
  # get the range  for the plot. If line two is bigger then line 1, the axis should be the size of line 2
  limx <- c(min(cc1$x, cc2$x, na.rm = T), max(cc1$x, cc2$x, na.rm = T))
  limy <- c(min(cc1$y, cc2$y, na.rm = T), max(cc1$y, cc2$y, na.rm = T))
  
  plot(cc1, type="b", col="red", lwd=3, xlim = limx, ylim = limy)
  if (!is.null(l2)) {
    lines(cc2, type="b", col="blue", lwd=3)
  }
  if (!is.null(l3)) {
    lines(getCoordinates(l3), type="b", col="green", lwd=3)
  }
  if (!is.null(l4)) {
    lines(getCoordinates(l4), type="b", col="purple", lwd=3)
  } 
}

plotGraph <- function(g) {
  plot(g, vertex.color=igraph::degree(g)+1)
  legend("topright", legend = sort(unique(igraph::degree(g))), col=sort(unique(igraph::degree(g)+1)), pch=19, title = "degree")
}

plotNumberOfPointsInRadius <- function(e, r) {
  
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      coordinatesOfLine <- getCoordinatesOfALine(e[[edge]])
      numberOfCoordinates <- nrow(coordinatesOfLine)
      for (point in 1:numberOfCoordinates) {
        p <- createPoint(coordinatesOfLine[point,])
        plot(p, add=T, pch=19, col=getNumberOfPointsInRadiusFromPoint(p, e, r))
      }
    }
  }
  return()
}

plotEdges <- function(e) {
  alreadyPlotted <- FALSE
  for (i in 1:length(e)) {
    edge <- e[i]
    if (is.na(e[i])) {
      next 
    } else {
      if (!alreadyPlotted) {
        plot(createPoint(getExtentOfEdges(edges)), cex=1, col="white")
        #plot(createPoint(as.data.frame(t(e[[i]]@bbox))), cex=1, col="white")
        alreadyPlotted <- TRUE
      }
      plot(e[[i]], add=T, type="b", col="black", lwd=2)
      plot(getStartPointOfLine(e[[i]]), add=T, col=getVertexDegree(getStartPointOfLine(e[[i]]), e), pch=19)
      text(getCoordinatesOfPoints(getStartPointOfLine(e[[i]])), labels=i, pos=1)
    }
  }
}

plotEdge <- function(x, i, color) {
  plot(x, type="b", add=T, col=color)
  #text(getCoordinatesOfPoints(getStartPointOfLine(x)), labels=??, pos=1)
}

plotEdgesBuffer <- function(x, b) {
  plot(gBuffer(x, width = b), add=T)
}

# program

loadFiles <- function() {
  listOfFiles=tk_choose.files(caption = "Select files", filters = matrix(c("topografix gpx", ".gpx"), ncol=2))
  edges <- lapply(X = listOfFiles, FUN = function(x) createLine(getCoordinatesOfPoints(readOGR(x, layer="track_points"))))
  
  return(edges)
}

reduceIntersections <- function(e) {
  
  # Find intersections with other edges and cut both edges into non intersecting parts, removing both original edges
  for (edge in 1:(length(e)-1)) {
    # edge <- edge +1
    if (is.na(e[edge])) {
      next
    } else {
      for (dge in (edge+1):length(e)) { 
        # dge <- dge +1
        # dge <- edge +1
        if (is.na(e[dge])) {
          next 
        } else {
          doIntersect <- gIntersects(e[[edge]], e[[dge]])
          if (doIntersect) { # They intersect
            dd <- gIntersection(e[[edge]], e[[dge]], byid = T)
            # plot(dd, add=T, pch=19, col="green", cex=2)
            # plot(e[[edge]], add=T, col="red", lwd=2)
            # plot(e[[dge]], add=T, col="blue", lwd=2)             
            intersections <- getCoordinatesOfPoints(dd) 
            
            # After two lines have been split they still intersect at the splitpoint, but i shouldn't split them.
            # i have written own touch function on start or end points
            doTouch <- touchPointOnStartOrEndNode(l1 = e[[edge]], l2 = e[[dge]], p = intersections)
            if (doTouch) {
              # points are on each other start or endpoints
              next
            } else {
              e <- splitLine(i = edge, ps = intersections, e = e)
              # plotEdges(e)
              e <- splitLine(i = dge, ps = intersections, e = e)
              # plotEdges(e)
              #e <- cleanEdgeList(e = e)
              # plotEdges(e)
              #               fileName <- "Documents/scripts/trailNetwork/tmp/reduce.csv"
              #               if(file.exists(fileName)) {
              #                 file.remove(fileName)
              #               }
              #               
              #               write.table(x = data.frame(edge = edge, dge = dge, n = 1), file =fileName, append = T, row.names = F, col.names = F)
              print(paste(edge, dge, sep=" - "))
              e <- reduceIntersections(e = e)
              return(e)
            }
          } else {
            next
          }
        } 
      }
    }
  }    
  return(e)
}

parallelReduction <- function(e) {
  # find two edges with same start and endpoint close together and 'averages' both in a new edge, removing both original edges
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    #for (edge in 1:(numberOfEdges-1)) {
    if (is.na(e[edge])) {
      next
    } else {
      ps1 <- getStartPointOfLine(l = e[[edge]])
      
      # Dit is volgens mij beter
      # wellicht wel langzamer, maar in elk geval overzichterlijker
      index <- getIndicesOfEqualPoints(p = ps1, pos = "start", e = e)
      
      pe1 <- getEndPointOfLine(l = e[[edge]])
      indexEnd <- getIndicesOfEqualPoints(p = pe1, pos = "end", e = e)
      
      index <- index[!(index[index %in% indexEnd] %in% edge)] # See if index of startpoints is in index of endpoints and remove own index
      if (length(index) == 1) {
        indexOfLineToAverage <- index
        newLine <- averageLines(e[[edge]], e[[indexOfLineToAverage]])
        e <- addEdge(newLine, e = e)
        e <- removeEdge(edge, e = e)
        e <- removeEdge(indexOfLineToAverage, e = e)
      }
      
      #       for (dge in (edge+1):numberOfEdges) {
      #         if (is.na(e[dge])) {
      #           next
      #         } else {
      #           ps2 <- getStartPointOfLine(e[[dge]])
      #           if (gEquals(ps1, ps2)) {
      #             pe1 <- getEndPointOfLine(e[[edge]])
      #             pe2 <- getEndPointOfLine(e[[dge]])
      #             if (gEquals(pe1, pe2)) {
      #               # we can average!
      #               newLine <- averageLines(e[[edge]], e[[dge]])
      #               e <- addEdge(newLine, e = e)
      #               e <- removeEdge(edge, e = e)
      #               e <- removeEdge(dge, e = e)
      #             }
      #           }
      #         }
      #       }
    }
  }
  return(e)
}

serialReduction <- function(e) {
  # Finds two conecting edges on a non incident vertex and combines them into one, removing both original edges
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next
    } else {
      endPoint <- getEndPointOfLine(e[[edge]])
      vertexDegree <- getVertexDegree(endPoint , e = e)
      if (vertexDegree == 2) {
        lineId <- getLineWithStartPoint(endPoint, e = e)
        if (!is.null(lineId)) {
          
          
          #        stop(paste("Found vertex. Member of", edge, "and", lineId, sep=" "))
          newLine <- combineLines(l1 = e[[edge]], l2 = e[[lineId]])
          e <- addEdge(newLine, e=e)
          e <- removeEdge(edge, e=e)
          e <- removeEdge(lineId, e=e)
          e <- serialReduction(e = e)
          break
        } else {
          # end point zijn gelijk
          # zoek lijn met gelijk endpoint
        }
      } else {
        # see if start point has degree of two
        startPoint <- getStartPointOfLine(e[[edge]])
        vertexDegree <- getVertexDegree(endPoint , e = e)
        if (vertexDegree == 2) {
          lineId <- getLineWithStartPoint(startPoint, e = e)
          
          #        stop(paste("Found vertex. Member of", edge, "and", lineId, sep=" "))
          newLine <- combineLines(l1 = e[[edge]], l2 = reverseLine(e[[lineId]]) )# reverse line to get correct ordering, by creating the newLine, line will be sorted again in correct direction
          e <- addEdge(newLine, e=e)
          e <- removeEdge(edge, e=e)
          e <- removeEdge(lineId, e=e)
          e <- serialReduction(e = e)
          break
        }
      }
    }
  }
  return(e)
}

edgeContraction <- function(tol, e) {
  # find an edge with a degree one start or end node/vertex and removes it when it is smaller than tol
  numberOfEdges <- length(e)
  for (edge in 1:numberOfEdges) {
    if (is.na(e[edge])) {
      next 
    } else {
      vertexDegreeOfStartPoint <- getVertexDegree(p = getStartPointOfLine(l = e[[edge]]), e = e)
      vertexDegreeOfEndPoint <- getVertexDegree(p = getEndPointOfLine(l = e[[edge]]), e = e)
      if (vertexDegreeOfStartPoint == 1 | vertexDegreeOfEndPoint == 1) {
        # It is an end or start line
        lengthOfLine <- SpatialLinesLengths(SL = e[[edge]], longlat = F) * 1000 # in meters
        if (lengthOfLine < tol) {
          e <- removeEdge(i = edge, e = e)
        }
      }
    }
  }
  return(e)
}

trailnetwork <- function() {
  edges <- list()
  edges <- addEdge(createLine(x = c(1, 4, 6, 8, 8.5, 9, 13, 15, 16, 17, 19), 
                              y = c(2, 2, 4, 1, 2.0, 5, 8 , 6 , 4 , 2 , 1), doSort=TRUE), edges)
  edges <- addEdge(createLine(x = c(1, 2, 6, 7, 7.5, 8, 10 , 11, 12, 12.5, 14, 16, 18), 
                              y = c(5, 4, 2, 3, 9.0, 7, 5 , 4 , 7 , 8 , 9 , 2 , 3), doSort=TRUE), edges)
  edges <- addEdge(createLine(x = c(10, 1), 
                              y = c(10, 1), doSort=TRUE), edges)
  plot(edges[[1]], type="b", col="green", xlim=c(0,20), ylim=c(0,10))#, xlim=c(2,7))
  plot(edges[[2]], add=T, type="b", col="green")
  plot(edges[[3]], add=T, type="b", col="green")
  edges <- reduceIntersections(e = edges)
  #  edges <- cleanEdgeList(e = edges)
  edges <- parallelReduction(e = edges)
  # edges <- cleanEdgeList(e = edges)
  edges <- serialReduction(e = edges)
  #edges <- cleanEdgeList(e = edges)
  edges <- edgeContraction(tol = 2.5, e = edges)
  edges <- cleanEdgeList(e = edges)
  plotEdges(edges)
  #return(edges)
}

#trailnetwork()