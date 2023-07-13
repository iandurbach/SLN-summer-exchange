library(dplyr)
library(raster)
library(secr)
library(sf)
library(ggplot2)
library(secrdesign)
library(ggspatial)
library(elevatr)

# coordinate reference system to use
my_crs <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"

# (1) Read in survey region

## Read in from shape file and reproject
survey_area <- st_read("data/survey_region/survey_region.shp")
survey_area <- survey_area %>% st_zm() %>% st_transform(my_crs) %>% st_geometry()

## Check area in km2
st_area(survey_area) / 1000000

# (2) Potential design 1: regular grid

## Set random origin for grid (note this might put grid outside survey area, if so redo)
set.seed(4321)
random_origin <- st_coordinates(sf::st_sample(survey_area, size = 1))

## Make 7x7 grid with 3km spacing
grid <- make.grid(nx = 7, ny = 7, spacing = 3000, detector = "count", originxy = random_origin)

# (3) Potential design 2: algorithmic (min[n,r]-type) design

## Set up grid of potential camera locations 
cams <- st_make_grid(survey_area, what = "corners", cellsize = c(2000,2000))

plot(survey_area)
plot(cams, add = TRUE)

## Only keep those that fall inside survey region
cams <- cams[survey_area] 

plot(survey_area)
plot(cams, add = TRUE)

## Drop any points < 2000 or above 4500 elevation
cams_any <- cams %>% elevatr::get_elev_point(prj = st_crs(survey_area)$proj4string, src = "aws", overwrite = TRUE) 
cams_hab <- cams_any %>% dplyr::filter(elevation > 2000, elevation < 4500) %>% st_geometry()

plot(survey_area)
plot(st_geometry(cams_any), add = TRUE)
plot(cams_hab, add = TRUE, col = "red")

## Final potential camera locations for design
cams <- cams_hab

## Make an integration mesh/mask, buffer 20km, spacing = 2km
mesh_poly <- st_union(st_buffer(cams, dist = 20000))
mesh_grid <- st_make_grid(mesh_poly, what = "corners", cellsize = c(2000, 2000))
mesh_grid <- mesh_grid[mesh_poly]

plot(mesh_poly, axes = TRUE)
plot(mesh_grid, add = TRUE)
plot(st_geometry(survey_area), add = TRUE, fill = NA)
plot(cams, add= TRUE, col = "red")

## Turn potential camera locations into an secr 'traps' object (needed by GAoptim function)
alltraps_df <- data.frame(x = st_coordinates(cams)[,1], y = st_coordinates(cams)[,2]) 
alltraps <- read.traps(data = alltraps_df, detector = "count")

## Turn mesh into an secr 'mask' object (needed by GAoptim function)
mesh_df <- data.frame(x = st_coordinates(mesh_grid)[,1], 
                      y = st_coordinates(mesh_grid)[,2]) 
mesh <- read.mask(data = mesh_df)

## Run optimization function to produce design (run 5 times here for speed, use 50-100 for real runs)
mnr <- GAoptim(mask = mesh, 
               alltraps = alltraps, 
               ntraps = 50, 
               detectpar = list(lambda0 = 0.25, sigma = 8000), 
               D = 0.0001, 
               criterion = 4, 
               seed = 1234, 
               noccasions = 1, detectfn = c("HHN"), 
               ngen = 5, verbose = TRUE)

# (4) Potential design 3: algorithmic (P2-type) design

p2 <- GAoptim(mask = mesh, 
              alltraps = alltraps, 
              ntraps = 50, 
              detectpar = list(lambda0 = 0.25, sigma = 8000), 
              D = 0.0001, 
              criterion = 6, 
              seed = 1234, 
              noccasions = 1, detectfn = c("HHN"), 
              ngen = 5, verbose = TRUE)

# (5) Potential design 3: lacework design

## Turn survey_area into a SpatialPolygonsDataFrame (needed by make.lacework to remove points outside survey area)
survey_area_sp <- as(survey_area[[1]], 'Spatial')
survey_area_sp <- as(survey_area_sp, 'SpatialPolygonsDataFrame')

set.seed(1234)
lw <- make.lacework(region = survey_area_sp, 
                          spacing = c(15000, 3000),  
                          rotate = 45, 
                          detector = "count", keep.design = TRUE)

plot(lw)

## Manually remove some traps so end up with roughly the same number of traps as other designs
lw <- lw[(lw$x > 300000) & (lw$y < 4283000), ]

plot(lw)

## (6) Post-processing to make nice plots

## Extract camera locations into a data.frame

mnr_traps <- as.data.frame(mnr$optimaltraps)
p2_traps <- as.data.frame(p2$optimaltraps)
grid_traps <- as.data.frame(grid)
lw_traps <- as.data.frame(lw)

## Turn into sf objects

cams_sf <- st_as_sf(cams, coords = c("x", "y"), crs = my_crs) %>% st_sf() 
mnr_traps_sf <- st_as_sf(mnr_traps, coords = c("x", "y"), crs = my_crs) %>% st_sf() 
p2_traps_sf <- st_as_sf(p2_traps, coords = c("x", "y"), crs = my_crs) %>% st_sf() 
grid_traps_sf <- st_as_sf(grid_traps, coords = c("x", "y"), crs = my_crs) %>% st_sf() 
lw_traps_sf <- st_as_sf(lw_traps, coords = c("x", "y"), crs = my_crs) 

# (7) Make 2x2km grid centred on selected camera locations

cellSize <- 2000

## Grid cells selected by regular grid design
grid_grids <- st_buffer(grid_traps_sf, dist = cellSize/2, endCapStyle="SQUARE", nQuadSegs = 1)

## Grid cells selected by min(n,r) design

### First need a 2x2km grid centred on all possible camera locations

cams_grid <- (st_bbox(cams) + cellSize/2*c(-1,-1,1,1)) %>%
  st_make_grid(cellsize=c(cellSize, cellSize)) %>% st_sf()

plot(survey_area)
plot(cams_grid, add= TRUE)
plot(cams, add= TRUE)

### Now extract those cells that contain a selected camera location
mnr_grids <- cams_grid[mnr_traps_sf, ]

## Grid cells selected by P2 design (same process)

p2_grids <- cams_grid[p2_traps_sf, ]

## Grid cells selected by lacework design
lw_grids <- st_buffer(lw_traps_sf, dist = cellSize/2, endCapStyle="SQUARE", nQuadSegs = 1)

lw_grids <- (st_bbox(lw_traps_sf) + cellSize/2*c(-1,-1,1,1)) %>% st_make_grid(cellsize=c(cellSize, cellSize)) %>% st_sf()
lw_grids <- lw_grids[lw_traps_sf, ]

#grid_grids <- (st_bbox(grid_traps_sf) + cellSize/2*c(-1,-1,1,1)) %>% st_make_grid(cellsize=c(cellSize, cellSize)) %>% st_sf()
#grid_grids <- grid_grids[grid_traps_sf, ]

# (8) Plots use ggspatial

esri_topo <- paste0('https://services.arcgisonline.com/arcgis/rest/services/', 'World_Topo_Map/MapServer/tile/${z}/${y}/${x}.jpeg')

g0 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 11) +
  layer_spatial(survey_area, colour = "red", fill = NA) + 
  annotation_scale(location = "tl") + 
  theme_minimal()

g0
# ggsave("images/survey_region.png", g0, width=8, height=6, dpi = 300)

g1 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 11) +
  layer_spatial(survey_area, colour = "red", fill = NA) + 
  layer_spatial(cams_sf, colour = "black", fill = "black", alpha = 0.4) + 
  layer_spatial(mnr_grids, colour = "blue") + 
  layer_spatial(mnr_traps_sf, colour = "blue", fill = "blue") + 
  annotation_scale(location = "tl") + 
  theme_minimal()

g1
# ggsave("images/mnr_survey.png", g1, width=8, height=6, dpi = 300)

g2 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 11) +
  layer_spatial(survey_area, colour = "red", fill = NA) + 
  layer_spatial(cams_sf, colour = "black", fill = "black", alpha = 0.4) + 
  layer_spatial(p2_grids, colour = "blue") + 
  layer_spatial(p2_traps_sf, colour = "blue", fill = "blue") + 
  annotation_scale(location = "tl") + 
  theme_minimal()

g2
# ggsave("images/p2_survey.png", g2, width=8, height=6, dpi = 300)

g3 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 11) +
  layer_spatial(survey_area, colour = "red", fill = NA) + 
  layer_spatial(lw_grids, colour = "blue") + 
  layer_spatial(lw_traps_sf, colour = "blue", fill = "blue") + 
  annotation_scale(location = "tl") + 
  theme_minimal()

g3
# ggsave("imageslw_survey.png", g3, width=8, height=6, dpi = 300)

g4 <- ggplot() +
  annotation_map_tile(type = esri_topo, zoom = 11) +
  layer_spatial(survey_area, colour = "red", fill = NA) + 
  layer_spatial(grid_grids, colour = "blue") + 
  layer_spatial(grid_traps_sf, colour = "blue", fill = "blue") + 
  annotation_scale(location = "tl") + 
  theme_minimal()

g4
# ggsave("images/grid_survey.png", g4, width=8, height=6, dpi = 300)

# (9) Compare designs 

## Assess expected sample sizes
Enrm(traps = mnr$optimaltraps, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
Enrm(traps = p2$optimaltraps, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
Enrm(traps = grid, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
Enrm(traps = lw, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")

## Assess expected CV
minnrRSE(traps = mnr$optimaltraps, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
minnrRSE(traps = p2$optimaltraps, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
minnrRSE(traps = grid, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")
minnrRSE(traps = lw, mask = mesh, detectpar = list(lambda0 = 0.25, sigma = 8000), D = 1/10000, noccasions = 1, detectfn = "HHN")

# save(mnr, p2, grid, lw, mesh, file = "output/all-designs.Rdata")
