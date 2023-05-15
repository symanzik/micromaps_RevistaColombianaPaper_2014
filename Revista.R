### R Code for Revista Colombiana Article
###
### Title: Linked Micromap Plots for South America -
###   General Design Considerations and Specific Adjustments
### Revista Colombiana de Estadística (2014), Vol. 37, No. 2, pp. 451-469.
###
### DOI: http://dx.doi.org/10.15446/rce.v37n2spe.47949
###
### Authors: Juergen Symanzik, XiaoTian Dai, Marc H. Weber,
###    Quinn Payton & Michael G. McManus
###
### May 26, 2014
###
### Tested with R 3.1.0 with all R packages updated to the
###   most recent version (such as micromap 1.8 and ggplot2 1.0.0)
###   prior to the final test run on May 26, 2014
###
### Revised on March 5, 2015, to run under R 3.1.2
###   with all R packages updated to the most recent
###   version (such as micromap 1.9.2)
###
### Revised on May 10, 2023, to run under R 4.3.0
###   with all R packages updated to the most recent
###   version (such as micromap 1.9.7)
###
### This R code is split into eight main parts:
###   Part 1: General Settings and R Packages
###   Part 2: Shapefile Modification for Brazil
###   Part 3: Shapefile Modification for other Countries from South America
###   Part 4: Figure 2 with all 12 Countries from South America
###   Part 5: Figure 1 of Different Boundaries for Brazil
###   Part 6: Linked Micromap Plots for Brazil (Figures 5 & 6)
###   Part 7: Linked Micromap Plots for Argentina (Figure 3)
###   Part 8: Linked Micromap Plots for Uruguay (Figure 4)
###
### This R code should be understood as a working version that was developed
### to create the results and figures in the associated journal article.
### Future testing and generalizations will follow, with the ultimate goal
### to make this functionality available in the micromap R package
### (or as an independent R package).


### Part 1: General Settings and R Packages

# setwd("")

# If TRUE, use local shapefiles and web pages; otherwise download from web.
# Note that most underlying web pages have changed since 2014.
# Therefore, we have to work with local shapefiles and web pages.

useLocalFiles <- TRUE

# Original figures in 2014 were produced in postscript (ps) format to meet journal
# requirements. However, there were some limitations, e.g., only the 26-letter
# alphabet could be used with no special characters. This version of the code
# allows to toggle between ps (option [1]) and pdf (option [2]) output.
# The pdf option is set as default.

outputFormat <- c("ps", "pdf")[2]


library(micromap)
library(ggplot2)
library(grid)
library(maptools)
# library(sp)
# library(RColorBrewer)

library(lattice)
library(rgeos)
library(XML)

# install additional fonts for postscript files on Windows
#
# library(extrafont)

# font_import()
# loadfonts(device = "win")


Greys5 <- brewer.pal(5, "Greys")[5:1]
Greys6 <- brewer.pal(6, "Greys")[6:1]


### Part 2: Shapefile Modification for Brazil

# Rescaling and Shifting Area Functions
# ========================================================

# Functions and pseudo-code for implementing rescaling small
# map areas and shifting them within the map. Numerous diagnostive
# plots and output is included to allow interested readers
# to better understand the functionality of this code
# and adapt the code to their own data.


# Example shapefile data for Brazil from maps library

if (useLocalFiles) {
  print(load("Shapefiles/BRA_adm1.RData"))
} else {
  con <- url("http://gadm.org/data/rda/BRA_adm1.RData")
  print(load(con))
  close(con)
}

head(gadm@data)
plot(gadm)


AreaPercent <- function(x) {
  tot_area <- sum(sapply(slot(x, "polygons"), slot, "area"))
  sapply(slot(x, "polygons"), slot, "area") / tot_area * 100
}

# We need to determine (empirically) for some maps which subareas are
# too small to be seen on a micromap, e.g., when the subarea is less
# than 1% or 1/2% of the entire area (this is just a guess). This could be
# a default argument to the RescaleArea function (but can also be adjusted by a user).

RescaleArea <- function(x, unit, rescale_factor) {
  b <- gBuffer(x[unit, ], width = rescale_factor, byid = TRUE)
  return(b)
}


# For each such area, shift it in such a way that the midpoint becomes (0, 0),
# multiply all coordinates with a common multiplier > 1 to enlarge the area
# proportionally (so that this area would be 1% or 1/2% of the entire area -
# or whatever else is our critical cutoff above). Then shift to the new location.

ShiftArea <- function(region, subregion, latlon, offset) {
  for (i in 1:length(region@polygons[[subregion]]@Polygons)) {
    cds <- slot(slot(
      slot(region, "polygons")[[subregion]],
      "Polygons"
    )[[i]], "coords")
    # get coords for vertices and label points
    l1 <- slot(slot(
      slot(region, "polygons")[[subregion]],
      "Polygons"
    )[[i]], "labpt")
    l2 <- slot(
      slot(region, "polygons")[[subregion]],
      "labpt"
    )
    # shift vertices coords
    cds[, latlon] <- cds[, latlon] + offset
    # shift lapt1 (centerpoint)
    newl1 <- l1[latlon] + offset
    # shift lapt2 (centerpoint)
    newl2 <- l2[latlon] + offset
    # put these shifted points back into the polygon object
    slot(slot(
      slot(region, "polygons")[[subregion]],
      "Polygons"
    )[[i]], "coords") <- cds
    slot(slot(
      slot(region, "polygons")[[subregion]],
      "Polygons"
    )[[i]], "labpt") <- newl1
    slot(
      slot(region, "polygons")[[subregion]],
      "labpt"
    ) <- newl2
  }
  return(region)
}


# See how much space there is in the corners of the rectangle specified by the
# bounding box and the convex hull that surrounds the main geographic area.
# If there is not enough space in the corners, then one could possibly enlarge
# the bounding box.
# If there is enough space in the "corners" of the square, shift there.
# If not, increase the plot area to [-1.1, 1] x [-1, 1] or to [-1, 1.1] x [-1, 1],
# i.e., add some space to the left or right (depending on whether the midpoint of
# this area is <0 or >= 0).
# Shift the enlarged area there.
# The same could be done for areas far away. It may be necessary to resize such
# areas as well. In case of Alaska, the resize function should shrink the original area.

as.SpatialPolygons.bbox <- function(x) {
  # Create unprojected bbox as spatial object
  bboxMat <- rbind(
    c(bbox(x)["x", "min"], bbox(x)["y", "min"]),
    c(bbox(x)["x", "min"], bbox(x)["y", "max"]),
    c(bbox(x)["x", "max"], bbox(x)["y", "max"]),
    c(bbox(x)["x", "max"], bbox(x)["y", "min"]),
    c(bbox(x)["x", "min"], bbox(x)["y", "min"])
  )
  # clockwise, 5 points to close it
  bboxSP <- SpatialPolygons(list(Polygons(list(Polygon(bboxMat)), "bbox")),
    proj4string = CRS(proj4string(x))
  )
  return(bboxSP)
}

fullArea <- function(x) {
  # somehow compare bbox to convex hull to find appropriate free corners of frame
  conhull <- gConvexHull(x, byid = FALSE, id = NULL)
  return(conhull)
}

# Try applying functions
# Hard-wired approach to getting rid of smallest areas based on percent
# Project first to appropriate projection for SA
# Need to get rid of islands and simplify polygons

brazil <- spTransform(gadm, CRS("+proj=laea +lat_0=0 +lon_0=-80"))
# plot it
plot(brazil)

# Determine which areas are going to be too small to see, and then either
# dissolve these areas into neighboring polygons or shift and enlarge
AreaPercent(brazil)

# plot areas less than 1% in a different color to examine
plot(brazil[c(2, 7, 8, 15, 19, 20, 26), ], add = TRUE, col = "Blue")

##############################

# solution for Distrito Federal
brazil2 <- brazil

## Change IDs
brazil2$ID_1[brazil2$ID_1 %in% c(436, 438)] <- "438"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(-7), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

newarea <- RescaleArea(brazil, 7, 50000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the map
fedist <- gIntersection(brazil2[8, ], newarea, byid = TRUE)

## Plot the output
plot(fedist, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(fedist) <- "436"
fedistdata <- brazil[7, ]@data
row.names(fedistdata) <- "436"
fedist <- SpatialPolygonsDataFrame(fedist, fedistdata)

# and cut the dissolved area with the new area
goias <- gDifference(brazil2[8, ], fedist)
row.names(goias) <- "438"
goiasdata <- brazil2[8, ]@data
row.names(goiasdata) <- "438"
goias <- SpatialPolygonsDataFrame(goias, goiasdata)
brazil2 <- brazil2[c(-8), ]
brazil3 <- rbind(brazil2, goias, fedist)
plot(brazil3)

################################

# solution for Rio de Janeiro
brazil2 <- brazil3

## Change IDs
brazil2$ID_1[brazil2$ID_1 %in% c(448, 442)] <- "442"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1:6, 27, 7, 26, 8:16, 18:25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

rio <- thinnedSpatialPoly(brazil[19, ], tolerance = 18000)
plot(brazil[19, ])
plot(rio, col = "red", add = TRUE)
length(slot(rio, "polygons"))
newarea <- RescaleArea(rio, 1, 50000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
rio <- gIntersection(brazil2[13, ], newarea, byid = TRUE)

## Plot the output
plot(rio, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(rio) <- "448"
riodata <- brazil[19, ]@data
row.names(riodata) <- "448"
rio <- SpatialPolygonsDataFrame(rio, riodata)

# and cut the dissolved area with the new area
minas <- gDifference(brazil2[13, ], rio)
row.names(minas) <- "442"
minasdata <- brazil2[13, ]@data
row.names(minasdata) <- "442"
minas <- SpatialPolygonsDataFrame(minas, minasdata)
brazil2 <- brazil2[c(-13), ]

brazil3 <- rbind(brazil2, rio, minas)
plot(brazil3)

################################

# solution for Esprito Santo
brazil2 <- brazil3

## Change IDs
brazil2$ID_1[brazil2$ID_1 %in% c(437, 442)] <- "442"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1:7, 9:12, 27, 13:17, 26, 18:25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

santo <- thinnedSpatialPoly(brazil[8, ], tolerance = 18000)
plot(brazil[8, ])
plot(santo, col = "red", add = TRUE)
length(slot(santo, "polygons"))
newarea <- RescaleArea(santo, 1, 50000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
santo <- gIntersection(brazil2[12, ], newarea, byid = TRUE)

## Plot the output
plot(santo, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(santo) <- "437"
santodata <- brazil[8, ]@data
row.names(santodata) <- "437"
santo <- SpatialPolygonsDataFrame(santo, santodata)

# and cut the dissolved area with the new area
minas <- gDifference(brazil2[12, ], santo)
row.names(minas) <- "442"
minasdata <- brazil2[12, ]@data
row.names(minasdata) <- "442"
minas <- SpatialPolygonsDataFrame(minas, minasdata)
brazil2 <- brazil2[c(-12), ]

brazil3 <- rbind(brazil2, santo, minas)
plot(brazil3)

################################

# solution for Sergipe
brazil2 <- brazil3

## Change IDs
brazil2@data
brazil2$ID_1[brazil2$ID_1 %in% c(455, 434)] <- "434"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1:7, 26, 8:11, 27, 12:23, 25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

sergipe <- thinnedSpatialPoly(brazil[26, ], tolerance = 18000)
plot(brazil[26, ])
plot(sergipe, col = "red", add = TRUE)
length(slot(sergipe, "polygons"))
newarea <- RescaleArea(sergipe, 1, 120000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
sergipe <- gIntersection(brazil2[5, ], newarea, byid = TRUE)

## Plot the output
plot(sergipe, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(sergipe) <- "455"
sergipedata <- brazil[26, ]@data
row.names(sergipedata) <- "455"
sergipe <- SpatialPolygonsDataFrame(sergipe, sergipedata)

# and cut the dissolved area with the new area
bahia <- gDifference(brazil2[5, ], sergipe)
row.names(bahia) <- "434"
bahiadata <- brazil2[5, ]@data
row.names(bahiadata) <- "434"
bahia <- SpatialPolygonsDataFrame(bahia, bahiadata)
brazil2 <- brazil2[c(-5), ]

brazil3 <- rbind(brazil2, sergipe, bahia)
plot(brazil3)

################################

# solution for Rio Grande do Norte
brazil2 <- brazil3

## Change IDs
brazil2@data
brazil2$ID_1[brazil2$ID_1 %in% c(449, 435)] <- "435"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1:4, 27, 5:18, 20:24, 26, 25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

rionorte <- thinnedSpatialPoly(brazil[20, ], tolerance = 18000)
plot(brazil[20, ])
plot(rionorte, col = "red", add = TRUE)
length(slot(rionorte, "polygons"))
newarea <- RescaleArea(rionorte, 1, 70000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
rionorte <- gIntersection(brazil2[6, ], newarea, byid = TRUE)

## Plot the output
plot(rionorte, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(rionorte) <- "449"
rionortedata <- brazil[20, ]@data
row.names(rionortedata) <- "449"
rionorte <- SpatialPolygonsDataFrame(rionorte, rionortedata)

# and cut the dissolved area with the new area
ceara <- gDifference(brazil2[6, ], rionorte)
row.names(ceara) <- "435"
cearadata <- brazil2[6, ]@data
row.names(cearadata) <- "435"
ceara <- SpatialPolygonsDataFrame(ceara, cearadata)
brazil2 <- brazil2[c(-6), ]

brazil3 <- rbind(brazil2, rionorte, ceara)
plot(brazil3)

################################

# solution for Alagoas
brazil2 <- brazil3

## Change IDs
brazil2@data
brazil2$ID_1[brazil2$ID_1 %in% c(431, 455)] <- "455"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1, 3:5, 27, 6:18, 26, 19:25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

alagoas <- thinnedSpatialPoly(brazil[2, ], tolerance = 18000)
plot(brazil[2, ])
plot(alagoas, col = "red", add = TRUE)
length(slot(alagoas, "polygons"))
newarea <- RescaleArea(alagoas, 1, 70000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
alagoas <- gIntersection(brazil2[25, ], newarea, byid = TRUE)

## Plot the output
plot(alagoas, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(alagoas) <- "431"
alagoasdata <- brazil[2, ]@data
row.names(alagoasdata) <- "431"
alagoas <- SpatialPolygonsDataFrame(alagoas, alagoasdata)

# and cut the dissolved area with the new area
sergipe <- gDifference(brazil2[25, ], alagoas)
row.names(sergipe) <- "455"
sergipedata <- brazil2[25, ]@data
row.names(sergipedata) <- "455"
sergipe <- SpatialPolygonsDataFrame(sergipe, sergipedata)
brazil2 <- brazil2[c(-25), ]

brazil3 <- rbind(brazil2, alagoas, sergipe)
plot(brazil3)

################################

# solution for Paraiba
brazil2 <- brazil3

## Change IDs
brazil2@data
brazil2$ID_1[brazil2$ID_1 %in% c(444, 446)] <- "446"

## Merge Polygons
brazil2.sp <- unionSpatialPolygons(brazil2, brazil2$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil2.sp, "polygons"), function(x) slot(x, "ID"))
brazil2@data[, 4]

## Merge data
brazil2.data <- brazil2@data[c(1, 26, 2:13, 15:24, 27, 25), ]
row.names(brazil2.data) <- row.names(brazil2.sp)

## Build the new SpatialPolygonsDataFrame
brazil2 <- SpatialPolygonsDataFrame(brazil2.sp, brazil2.data)
plot(brazil2)

paraiba <- thinnedSpatialPoly(brazil[15, ], tolerance = 20000)
plot(brazil[15, ])
plot(paraiba, col = "red", add = TRUE)
length(slot(paraiba, "polygons"))
newarea <- RescaleArea(paraiba, 1, 40000)
plot(newarea, col = "red", add = TRUE)

# add enlarged area back in
## Clip the new area
paraiba <- gIntersection(brazil2[16, ], newarea, byid = TRUE)

## Plot the output
plot(paraiba, col = "khaki", bg = "azure2", add = TRUE)

## looks good, add back in
row.names(paraiba) <- "444"
paraibadata <- brazil[15, ]@data
row.names(paraibadata) <- "444"
paraiba <- SpatialPolygonsDataFrame(paraiba, paraibadata)

# and cut the dissolved area with the new area
pernambuco <- gDifference(brazil2[16, ], paraiba)
row.names(pernambuco) <- "446"
pernambucodata <- brazil2[16, ]@data
row.names(pernambucodata) <- "446"
pernambuco <- SpatialPolygonsDataFrame(pernambuco, pernambucodata)
brazil2 <- brazil2[c(-16), ]

brazil3 <- rbind(brazil2, paraiba, pernambuco)
plot(brazil3)

brazil4 <- thinnedSpatialPoly(brazil3, tolerance = 18000, minarea = 500)
plot(brazil4)

# try pushing federal district out into ocean to demonstrate shifting
brazil5 <- ShiftArea(brazil4, 7, 1, 1000000)
brazil5 <- ShiftArea(brazil5, 7, 2, -600000)
plot(brazil5)
fedist <- brazil5[7, ]

## Change IDs
brazil3$ID_1[brazil3$ID_1 %in% c(436, 438)] <- "438"

## Merge Polygons
brazil3.sp <- unionSpatialPolygons(brazil3, brazil3$ID_1)

## Rownames of the associated data frame must be the same as polygons IDs
sapply(slot(brazil3.sp, "polygons"), function(x) slot(x, "ID"))
brazil3@data[, 4]

## Merge data
brazil3.data <- brazil3@data[c(1, 26, 2:5, 7:24, 27, 25), ]
row.names(brazil3.data) <- row.names(brazil3.sp)

## Build the new SpatialPolygonsDataFrame
brazil3 <- SpatialPolygonsDataFrame(brazil3.sp, brazil3.data)
plot(brazil3)
brazil3 <- rbind(brazil3, fedist)
brazil3 <- thinnedSpatialPoly(brazil3, tolerance = 18000, minarea = 500)
plot(brazil3)


### Part 3: Shapefile Modification for other Countries from South America

ArgShapefile <- readShapeSpatial("Shapefiles/ARG_adm/ARG_adm1",
  verbose = TRUE
)
plot(ArgShapefile)
ArgShapefileThin <- thinnedSpatialPoly(ArgShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(ArgShapefileThin)

BolShapefile <- readShapeSpatial("Shapefiles/BOL_adm/BOL_adm1",
  verbose = TRUE
)
plot(BolShapefile)
BolShapefileThin <- thinnedSpatialPoly(BolShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(BolShapefileThin)

BraShapefile <- readShapeSpatial("Shapefiles/BRA_adm/BRA_adm1",
  verbose = TRUE
)
plot(BraShapefile)
BraShapefileThin <- thinnedSpatialPoly(BraShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(BraShapefileThin)

ChiShapefile <- readShapeSpatial("Shapefiles/CHL_adm/CHL_adm1",
  verbose = TRUE
)
plot(ChiShapefile)
ChiShapefileThin <- thinnedSpatialPoly(ChiShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(ChiShapefileThin, xlim = c(-90, -50))

ColShapefile <- readShapeSpatial("Shapefiles/COL_adm/COL_adm1",
  verbose = TRUE
)
plot(ColShapefile)
ColShapefileThin <- thinnedSpatialPoly(ColShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(ColShapefileThin)

EcuShapefile <- readShapeSpatial("Shapefiles/ECU_adm/ECU_adm1",
  verbose = TRUE
)
plot(EcuShapefile)
EcuShapefileThin <- thinnedSpatialPoly(EcuShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(EcuShapefileThin)

GuyShapefile <- readShapeSpatial("Shapefiles/GUY_adm/GUY_adm1",
  verbose = TRUE
)
plot(GuyShapefile)
GuyShapefileThin <- thinnedSpatialPoly(GuyShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(GuyShapefileThin)

PerShapefile <- readShapeSpatial("Shapefiles/PER_adm/PER_adm1",
  verbose = TRUE
)
plot(PerShapefile)
PerShapefileThin <- thinnedSpatialPoly(PerShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(PerShapefileThin)

PryShapefile <- readShapeSpatial("Shapefiles/PRY_adm/PRY_adm1",
  verbose = TRUE
)
plot(PryShapefile)
PryShapefileThin <- thinnedSpatialPoly(PryShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(PryShapefileThin)

SurShapefile <- readShapeSpatial("Shapefiles/SUR_adm/SUR_adm1",
  verbose = TRUE
)
plot(SurShapefile)
SurShapefileThin <- thinnedSpatialPoly(SurShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(SurShapefileThin)

UryShapefile <- readShapeSpatial("Shapefiles/URY_adm/URY_adm1",
  verbose = TRUE
)
plot(UryShapefile)
UryShapefileThin <- thinnedSpatialPoly(UryShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(UryShapefileThin)

VenShapefile <- readShapeSpatial("Shapefiles/VEN_adm/VEN_adm1",
  verbose = TRUE
)
plot(VenShapefile)
VenShapefileThin <- thinnedSpatialPoly(VenShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)
plot(VenShapefileThin)


### Part 4: Figure 2 with all 12 Countries from South America

if (outputFormat == "ps") {
  postscript("Output/SA_Countries_05_10_2023.ps",
    paper = "special",
    height = 12, width = 7.5, horizontal = FALSE
  )
} else {
  pdf("Output/SA_Countries_05_10_2023.pdf",
    paper = "special",
    height = 12, width = 7.5
  )
}

par(mfrow = c(4, 3), mar = c(1, 0, 1, 0))

plot(ArgShapefileThin)
title(main = "Argentina")
plot(BolShapefileThin)
title(main = "Bolivia")
plot(BraShapefileThin)
title(main = "Brazil")
plot(ChiShapefileThin, xlim = c(-90, -50))
title(main = "Chile")
plot(ColShapefileThin)
title(main = "Colombia")
plot(EcuShapefileThin)
title(main = "Ecuador")
plot(GuyShapefileThin)
title(main = "Guyana")
plot(PryShapefileThin)
title(main = "Paraguay")
plot(PerShapefileThin)
title(main = "Peru")
plot(SurShapefileThin)
title(main = "Suriname")
plot(UryShapefileThin)
title(main = "Uruguay")
plot(VenShapefileThin)
title(main = "Venezuela")

dev.off()


### Part 5: Figure 1 of Different Boundaries for Brazil

BraShapefile <- readShapeSpatial("Shapefiles/BRA_adm/BRA_adm1",
  verbose = TRUE
)
BraShapefileThin <- thinnedSpatialPoly(BraShapefile,
  tolerance = 0.05,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)

if (outputFormat == "ps") {
  postscript("Output/Brazil_Boundaries_05_10_2023.ps",
    paper = "special",
    height = 7.5, width = 7.5, horizontal = FALSE
  )
} else {
  pdf("Output/Brazil_Boundaries_05_10_2023.pdf",
    paper = "special",
    height = 7.5, width = 7.5
  )
}

par(mfrow = c(2, 2), mar = c(0, 0, 0, 0))

plot(BraShapefile)
plot(BraShapefileThin)
plot(brazil4)
plot(brazil3)

dev.off()


### Part 6: Linked Micromap Plots for Brazil (Figures 5 & 6)

if (useLocalFiles) {
  Brafile <- "WebPages/States of Brazil - Wikipedia, the free encyclopedia.htm"
  Bratable <- readHTMLTable(Brafile)[[12]]
} else {
  Braurl <- "http://en.wikipedia.org/wiki/States_of_Brazil"
  Bratable <- readHTMLTable(Braurl)[[12]]
}

Bratable <- apply(Bratable, 2, function(x) gsub("[(].*", "", x))
Bratable <- apply(Bratable, 2, function(x) gsub(".*00000", "", x))
Bratable <- apply(Bratable, 2, function(x) gsub("[,]", "", x))
Bratable <- apply(Bratable, 2, function(x) gsub("[%]", "", x))
Bratable[, 12] <- substr(Bratable[, 12], 1, nchar(Bratable[, 12]) - 1)
Bratable <- as.data.frame(Bratable)
Bratable <- data.frame(
  Bratable[, 1:4],
  sapply(Bratable[, 5:13], function(x) as.numeric(as.character(x)))
)
State <- as.character(sort(Bratable[, 2]))

###

StateReordered <- State[c(16, 18:27, 15, 17)]
State[15:27] <- StateReordered

brazil4@data <- cbind(brazil4@data, State)
Bratable[, 6] <- Bratable[, 6] / 1000000

BrazilPolys <- create_map_table(brazil4, "State")

Bratable$StateMod <- as.character(Bratable$State)

### .ps format
if (outputFormat == "ps") {
  Bratable$StateMod[3] <- "Amapa"
  Bratable$StateMod[6] <- "Ceara"
  Bratable$StateMod[8] <- "Espirito Santo"
  Bratable$StateMod[9] <- "Goias"
  Bratable$StateMod[10] <- "Maranhao"
  Bratable$StateMod[14] <- "Para"
  Bratable$StateMod[15] <- "Paraiba"
  Bratable$StateMod[16] <- "Parana"
  Bratable$StateMod[18] <- "Piaui"
  Bratable$StateMod[22] <- "Rondonia"
  Bratable$StateMod[25] <- "Sao Paulo"
}


# bra1_1 - figure not used in accompanying article
mmplot(
  stat.data = Bratable,
  map.data = BrazilPolys,
  panel.types = c(
    "map", "dot_legend", "labels",
    "dot", "dot", "dot"
  ),
  panel.data = list(
    NA, "points", "StateMod",
    "Population..2014.", "GDP.per.capita..R....2014.",
    "Infant.mortality..2014."
  ),
  ord.by = "GDP.per.capita..R....2014.",
  rev.ord = TRUE,
  grouping = c(4, 4, 4, 3, 4, 4, 4),
  vertical.align = "center",
  median.row = FALSE,
  map.link = c("State", "ID"),
  plot.width = 8,
  plot.panel.spacing = 1.5,
  colors = Greys5,
  map.color2 = gray(0.9),
  panel.att = list(
    list(1,
      header = "Maps",
      panel.width = 0.6,
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5
    ),
    list(2, point.size = 0.8),
    list(3,
      header = "States",
      align = "left",
      text.size = 0.9,
      panel.width = 1.3
    ),
    list(4,
      header = "Population\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 20, 40),
      xaxis.labels = list(0, 20, 40),
      xaxis.title = "People (Millions)"
    ),
    list(5,
      header = "GDP per Capita\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 25000, 50000),
      xaxis.labels = list(0, 25000, 50000),
      xaxis.title = "Brazillian Real ($)"
    ),
    list(6,
      header = "Infant Mortality\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 10, 20),
      xaxis.labels = list(0, 10, 20),
      xaxis.title = "Per thousand (per mille)"
    )
  ),
  print.file = paste0("Output/BRA_1_1_05_10_2023_not_used.", outputFormat),
  print.res = 300
)

# bra1_2 - Figure 5
mmplot(
  stat.data = Bratable,
  map.data = BrazilPolys,
  panel.types = c(
    "dot_legend", "labels", "dot",
    "dot", "dot", "map"
  ),
  panel.data = list(
    "points", "StateMod",
    "Population..2014.", "GDP.per.capita..R....2014.",
    "Infant.mortality..2014.", NA
  ),
  ord.by = "GDP.per.capita..R....2014.",
  rev.ord = TRUE,
  grouping = c(4, 4, 4, 3, 4, 4, 4),
  vertical.align = "center",
  median.row = FALSE,
  map.link = c("State", "ID"),
  plot.width = 8,
  plot.panel.spacing = 1.5,
  colors = Greys5,
  map.color2 = gray(0.9),
  panel.att = list(
    list(6,
      header = "Maps",
      panel.width = 0.6,
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5
    ),
    list(1, point.size = 0.8),
    list(2,
      header = "States",
      align = "left",
      text.size = 0.9,
      panel.width = 1.0
    ),
    list(3,
      header = "Population\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 20, 40),
      xaxis.labels = list(0, 20, 40),
      xaxis.title = "People (Millions)"
    ),
    list(4,
      header = "GDP per Capita\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 25000, 50000),
      xaxis.labels = list(0, 25000, 50000),
      xaxis.title = "Brazillian Real ($)"
    ),
    list(5,
      header = "Infant Mortality\nin 2014",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 10, 20),
      xaxis.labels = list(0, 10, 20),
      xaxis.title = "Per thousand (per mille)"
    )
  ),
  print.file = paste0("Output/BRA_1_2_05_10_2023.", outputFormat),
  print.res = 300
)

# bra2 - Figure 6

if (useLocalFiles) {
  Brafile2 <- "WebPages/List of Brazilian states by murder rate - Wikipedia, the free encyclopedia.htm"
  Bratable2 <- readHTMLTable(Brafile2)[[1]]
} else {
  Braurl2 <- "http://en.wikipedia.org/wiki/List_of_Brazilian_states_by_murder_rate"
  Bratable2 <- readHTMLTable(Braurl2)[[1]]
}

Bratable2 <- Bratable2[-c(8, 18, 23, 27, 32, 33), ]
Bratable2[, 2:13] <- sapply(Bratable2[, 2:13], function(x) as.numeric(as.character(x)))
Bratable2[, 1] <- as.character(Bratable2[, 1])
Bratable2[, 1] <- substring(Bratable2[, 1], 2)
Bratable2[, "change"] <- Bratable2[, 12] - Bratable2[, 2]
Bratable2[, "zero"] <- 0
Bratable2 <- Bratable2[order(Bratable2[, 1]), ]
myTable <- create_DF_rank(Bratable2,
  ord.by = "2008",
  group = c(4, 4, 4, 3, 4, 4, 4), rev.ord = TRUE
)
myTable[13:15, "pGrpOrd"] <- myTable[13:15, "pGrpOrd"] + 0.6
names(Bratable2)[1] <- "states"
names(myTable)[2] <- "states"

myAtts <- sample_att()
myNumber <- 1
myAtts$colors <- c("red", "orange", "green", "blue")
myAtts[[myNumber]]$panel.data <- c("zero", "change")
myColors <- myAtts$colors
myColumns <- myAtts[[myNumber]]$panel.data
myTable$data1 <- myTable[, myColumns[1]]
myTable$data2 <- myTable[, myColumns[2]]
myPanel <- ggplot(myTable) +
  geom_segment(
    aes(
      x = data1, y = -pGrpOrd,
      xend = data2, yend = -pGrpOrd, colour = factor(color)
    ),
    arrow = arrow(length = unit(0.1, "cm"))
  ) +
  facet_grid(pGrp ~ .) +
  scale_colour_manual(values = myColors, guide = "none")

assimilatePlot(myPanel, myNumber, myAtts)
myPanelAtts <- standard_att()
myPanelAtts <- append(
  myPanelAtts,
  list(line.width = 1, tip.length = 1)
)

arrow_plot_att <- function() {
  myPanelAtts <- standard_att()
  myPanelAtts <- append(
    myPanelAtts,
    list(line.width = 1, tip.length = 1)
  )
}

arrow_plot_build <- function(myPanel, myNumber, myNewStats, myAtts) {
  myColors <- myAtts$colors
  myColumns <- myAtts[[myNumber]]$panel.data
  myLineWidth <- myAtts[[myNumber]]$line.width
  myTipLength <- myAtts[[myNumber]]$tip.length
  myTable$data1 <- myNewStats[, myColumns[1]]
  myTable$data2 <- myNewStats[, myColumns[2]]
  myPanel <- ggplot(myTable) +
    geom_segment(
      aes(
        x = data1, y = -pGrpOrd,
        xend = data2, yend = -pGrpOrd,
        colour = factor(color)
      ),
      arrow = arrow(length = unit(0.1 * myTipLength, "cm")),
      size = myLineWidth
    ) +
    facet_grid(pGrp ~ ., space = "free", scales = "free_y") +
    scale_colour_manual(values = myColors, guide = "none")
  myPanel <- assimilatePlot(myPanel, myNumber, myAtts)
}

Bratable2$StateMod <- as.character(Bratable2$states)

### .ps format
if (outputFormat == "ps") {
  Bratable2$StateMod[3] <- "Amapa"
  Bratable2$StateMod[6] <- "Ceara"
  Bratable2$StateMod[8] <- "Espirito Santo"
  Bratable2$StateMod[9] <- "Goias"
  Bratable2$StateMod[10] <- "Maranhao"
  Bratable2$StateMod[14] <- "Para"
  Bratable2$StateMod[15] <- "Paraiba"
  Bratable2$StateMod[16] <- "Parana"
  Bratable2$StateMod[18] <- "Piaui"
  Bratable2$StateMod[22] <- "Rondonia"
  Bratable2$StateMod[25] <- "Sao Paulo"
}

###

mmplot(
  stat.data = Bratable2,
  map.data = BrazilPolys,
  panel.types = c("map", "dot_legend", "labels", "dot", "arrow_plot"),
  panel.data = list(NA, "points", "StateMod", "1998", c("zero", "change")),
  ord.by = "1998",
  rev.ord = TRUE,
  grouping = c(4, 4, 4, 3, 4, 4, 4),
  vertical.align = "center",
  map.link = c("states", "ID"),
  colors = Greys5,
  map.color2 = gray(0.9),
  panel.att = list(
    list(1,
      header = "Maps", panel.width = 0.8,
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5
    ),
    list(2, point.size = 0.8),
    list(3,
      header = "States",
      align = "left",
      text.size = 0.9,
      panel.width = 1.3
    ),
    list(4,
      header = "Murder Rate\nin 1998",
      point.type = 20,
      point.size = 1.4,
      graph.bgcolor = gray(0.9),
      xaxis.ticks = list(0, 20, 40, 60),
      xaxis.labels = list(0, 20, 40, 60),
      xaxis.title = "Percent (%)",
      panel.width = 1.8
    ),
    list(5,
      header = "Murder Rate Change\nfrom 1998 to 2008",
      graph.bgcolor = gray(0.9),
      xaxis.title = "Percentage Points",
      panel.width = 1.8
    )
  ),
  print.file = paste0("Output/BRA_2_05_10_2023.", outputFormat),
  print.res = 300
)


### Part 7: Linked Micromap Plots for Argentina (Figure 3)

# read in shape file
ArgShapefile <- readShapeSpatial("Shapefiles/ARG_adm/ARG_adm1",
  verbose = TRUE
)
names(ArgShapefile)
ArgShapefile$NAME_1
ArgShapefileThin <- thinnedSpatialPoly(ArgShapefile,
  tolerance = 0.1,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)

# manipulate population dataset
if (useLocalFiles) {
  Argfile <- "WebPages/List of Argentine provinces by population - Wikipedia, the free encyclopedia.htm"
  Argtable <- readHTMLTable(Argfile)[[1]]
} else {
  Argurl <- "http://en.wikipedia.org/wiki/List_of_Argentine_provinces_by_population"
  Argtable <- readHTMLTable(Argurl)[[1]]
}

Argtable <- Argtable[, -3]
Argtable <- apply(Argtable, 2, function(x) gsub(".*00000", "", x))
Argtable <- apply(Argtable, 2, function(x) gsub(",", "", x))
Argtable <- as.data.frame(Argtable)
Argtable[, 1] <- as.character(Argtable[, 1])
Argtable[, 2] <- as.character(Argtable[, 2])
Argtable[, 3] <- as.numeric(as.character(Argtable[, 3]))
Argtable[, 4] <- as.numeric(as.character(Argtable[, 4]))
Argtable[4, 2] <- "Ciudad de Buenos Aires"
provincematch <- gregexpr("\\<[[:alpha:][:space:]]*", Argtable[, 2])
Argtable[, 2] <- unlist(regmatches(Argtable[, 2], provincematch))
colnames(Argtable) <- c("Rank", "provinces", "2010 Population", "2001 Population")
Argtable[, "change"] <- Argtable[, 3] - Argtable[, 4]
provinces <- as.character(sort(Argtable[, 2]))
ArgShapefileThin@data <- cbind(ArgShapefileThin@data, provinces)
Argtable[, 3:5] <- sapply(Argtable[, 3:5], function(x) x / 1000000)
Argtable[, "ratio"] <- Argtable[, 3] / Argtable[, 4]

# arg1_1 - figure not used in accompanying article

ArgPolys <- create_map_table(ArgShapefileThin, "provinces")

myTable <- create_DF_rank(Argtable,
  ord.by = "2010 Population",
  group = c(5, 5, 4, 5, 5), rev.ord = TRUE
)
myTable[11:14, "pGrpOrd"] <- myTable[11:14, "pGrpOrd"] + 0.5
myAtts <- sample_att()
myNumber <- 1
myAtts$colors <- c("red", "orange", "green", "blue", "purple")
myAtts[[myNumber]]$panel.data <- "change"
myColors <- myAtts$colors
myColumns <- myAtts[[myNumber]]$panel.data
myTable$data1 <- myTable[, myColumns]
myPanel <- ggplot(myTable) +
  geom_segment(aes(
    x = 0, y = -pGrpOrd, xend = data1, yend = -pGrpOrd,
    colour = factor(color)
  )) +
  facet_grid(pGrp ~ .) +
  scale_colour_manual(
    values = myColors,
    guide = "none"
  )
assimilatePlot(myPanel, myNumber, myAtts)

myPanelAtts <- standard_att()
myPanelAtts <- append(myPanelAtts, list(line.width = 4))

bar_plot_att <- function() {
  myPanelAtts <- standard_att()
  myPanelAtts <- append(myPanelAtts, list(line.width = 4))
}

bar_plot_build <- function(myPanel, myNumber, myNewStats, myAtts) {
  myColors <- myAtts$colors
  myColumns <- myAtts[[myNumber]]$panel.data
  myLineWidth <- myAtts[[myNumber]]$line.width
  myTipLength <- myAtts[[myNumber]]$tip.length
  myTable$data1 <- myNewStats[, myColumns]
  myPanel <- ggplot(myTable) +
    geom_segment(
      aes(
        x = 0, y = -pGrpOrd,
        xend = data1, yend = -pGrpOrd,
        colour = factor(color)
      ),
      size = myLineWidth
    ) +
    facet_grid(pGrp ~ ., space = "free", scales = "free_y") +
    scale_colour_manual(values = myColors, guide = "none")
  myPanel <- assimilatePlot(myPanel, myNumber, myAtts)
}

Argtable$ProvinceTemp <- as.character(Argtable$provinces)

### .ps format
if (outputFormat == "ps") {
  Argtable$ProvinceTemp[2] <- "Cordoba"
  Argtable$ProvinceTemp[6] <- "Tucuman"
  Argtable$ProvinceTemp[8] <- "Entre Rios"
  Argtable$ProvinceTemp[15] <- "Rio Negro"
  Argtable$ProvinceTemp[16] <- "Neuquen"
}

###

mmplot(
  stat.data = Argtable,
  map.data = ArgPolys,
  panel.types = c("map", "labels", "dot", "bar_plot"),
  panel.data = list(NA, "ProvinceTemp", "2010 Population", "change"),
  ord.by = "2010 Population",
  rev.ord = TRUE,
  grouping = c(5, 5, 4, 5, 5),
  vertical.align = "center",
  map.link = c("provinces", "ID"),
  colors = Greys6,
  map.color2 = gray(0.95),
  panel.att = list(
    list(1,
      header = "Maps",
      panel.width = 0.6,
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5
    ),
    list(2,
      header = "Provinces",
      align = "center",
      text.size = 0.9,
      panel.width = 1.5
    ),
    list(3,
      header = "2010 Population",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(0, 5, 10, 15, 20),
      xaxis.labels = list(0, 5, 10, 15, 20),
      xaxis.title = "People (Millions)",
      panel.width = 1.6
    ),
    list(4,
      header = "Population Increase\nfrom 2001 to 2010",
      graph.bgcolor = gray(0.95),
      xaxis.title = "People (Millions)",
      panel.width = 1.6
    )
  ),
  print.file = paste0("Output/ARG_1_1_05_10_2023_not_used.", outputFormat),
  print.res = 300
)

# arg1_2 - Figure 3
mmplot(
  stat.data = Argtable,
  map.data = ArgPolys,
  panel.types = c(
    "map", "dot_legend", "labels",
    "dot", "bar_plot", "dot"
  ),
  panel.data = list(
    NA, "points", "ProvinceTemp",
    "2010 Population", "change", "ratio"
  ),
  ord.by = "2010 Population",
  rev.ord = TRUE,
  grouping = c(5, 5, 4, 5, 5),
  vertical.align = "center",
  map.link = c("provinces", "ID"),
  colors = Greys6,
  map.color2 = gray(0.95),
  plot.width = 10,
  panel.att = list(
    list(1,
      header = "Maps",
      panel.width = 0.6,
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5
    ),
    list(2, point.size = 0.8),
    list(3,
      header = "Provinces",
      align = "left",
      text.size = 1.0,
      panel.width = 1.5
    ),
    list(4,
      header = "2010 Population",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(0, 4, 8, 12, 16),
      xaxis.labels = list(0, 4, 8, 12, 16),
      xaxis.title = "People (Millions)",
      panel.width = 1.8
    ),
    list(5,
      header = "Population Increase\nfrom 2001 to 2010",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(0, 0.9, 1.8),
      xaxis.labels = list(0, 0.9, 1.8),
      xaxis.title = "People (Millions)",
      panel.width = 1.8
    ),
    list(6,
      header = "Ratio of 2010 Population\nto 2001 Population",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(1, 1.2, 1.4),
      xaxis.labels = list(1, 1.2, 1.4),
      panel.width = 1.8
    )
  ),
  print.file = paste0("Output/ARG_1_2_05_10_2023.", outputFormat),
  print.res = 300
)


### Part 8: Linked Micromap Plots for Uruguay (Figure 4)

# read in shape file
UruShapefile <- readShapeSpatial("Shapefiles/URY_adm/URY_adm1",
  verbose = TRUE
)
names(UruShapefile)
UruShapefile$NAME_1
UruShapefileThin <- thinnedSpatialPoly(UruShapefile,
  tolerance = 0.01,
  topologyPreserve = TRUE, avoidGEOS = FALSE
)

# manipulate population dataset

if (useLocalFiles) {
  Urufile <- "WebPages/Uruguay - Wikipedia, the free encyclopedia.htm"
  Urutable <- readHTMLTable(Urufile)[[3]]
} else {
  Uruurl <- "http://en.wikipedia.org/wiki/Uruguay"
  Urutable <- readHTMLTable(Uruurl)[[3]]
}

Urutable <- Urutable[-c(1, 21:22), ]
Urutable <- apply(Urutable, 2, function(x) gsub(".*00000", "", x))
Urutable <- apply(Urutable, 2, function(x) gsub(",", "", x))

Urutable <- data.frame(Urutable)
Urutable[, 1] <- as.character(Urutable[, 1])
Urutable[, 3] <- as.numeric(as.character(Urutable[, 3]))
Urutable[, 5] <- as.numeric(as.character(Urutable[, 5]))
names(Urutable)[1] <- "Department"
names(Urutable)[3] <- "Area"
names(Urutable)[5] <- "Population"

department <- sort(Urutable[, 1])
UruShapefileThin@data <- cbind(UruShapefileThin@data, department)
Urutable[, 3] <- Urutable[, 3] / 1000
Urutable[, 5] <- Urutable[, 5] / 1000

# uru1 - Figure 4

UruPolys <- create_map_table(UruShapefileThin, "department")

Urutable$deptTemp <- as.character(Urutable$Department)

### .ps format
if (outputFormat == "ps") {
  Urutable$deptTemp[11] <- "Paysandu"
  Urutable$deptTemp[12] <- "Rio Negro"
  Urutable$deptTemp[16] <- "San Jose"
  Urutable$deptTemp[18] <- "Tacuarembo"
}

###

mmplot(
  stat.data = Urutable,
  map.data = UruPolys,
  panel.types = c("dot_legend", "labels", "dot", "dot", "map"),
  panel.data = list("labels", "deptTemp", "Area", "Population", NA),
  ord.by = "Area",
  rev.ord = TRUE,
  grouping = c(5, 4, 4, 5),
  vertical.align = "center",
  median.row = TRUE,
  two.ended = TRUE,
  map.link = c("Department", "ID"),
  plot.width = 6,
  colors = Greys6,
  map.color2 = gray(0.95),
  panel.att = list(
    list(5,
      header = "Maps",
      active.border.color = "black",
      active.border.size = 0.5,
      inactive.border.color = gray(0.7),
      inactive.border.size = 0.5,
      panel.width = 1.8
    ),
    list(1, point.size = 0.8),
    list(2,
      header = "Department",
      align = "left",
      text.size = 1.0,
      panel.width = 1.2
    ),
    list(3,
      header = "Area",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(0, 10, 20),
      xaxis.labels = list(0, 10, 20),
      xaxis.title = "km2 (Thousands)",
      panel.width = 1.4
    ),
    list(4,
      header = "2011 Population",
      graph.bgcolor = gray(0.95),
      xaxis.ticks = list(0, 650, 1300),
      xaxis.labels = list(0, 650, 1300),
      xaxis.title = "People (Thousands)",
      panel.width = 1.4
    )
  ),
  print.file = paste0("Output/URU_1_05_10_2023.", outputFormat),
  print.res = 300
)
