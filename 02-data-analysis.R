library(dplyr)
library(raster)
library(secr)
library(sf)
library(ggplot2)
library(secrdesign)
library(ggspatial)

#load("output/design_outputs.Rdata")

my_crs <- "+proj=utm +zone=43 +datum=WGS84 +units=m +no_defs"

# (0) Preliminary stuff

## Read in survey region again
survey_area <- st_read("data/survey_region/survey_region.shp")
survey_area <- survey_area %>% st_zm() %>% st_transform(my_crs) %>% st_geometry()

# (1) Read in trapfile

## Read in
traps_df <- read.csv("data/secr_inputs/ex_trapfile_ll.csv") %>% as.data.frame()

## Convert lat long to UTM
traps_sf <- st_as_sf(traps_df, coords = c("x", "y")) %>% st_set_crs(4326) %>% st_transform(my_crs)
traps_df$x <- st_coordinates(traps_sf)[,1]
traps_df$y = st_coordinates(traps_sf)[,2]

## Only keep trapID, x, y, effort/usage 
## This is because you can't add covariates directly if using read.traps with "data" (frame) option
traps_df_nocovs <- traps_df[,1:4]
traps <- read.traps(data = traps_df_nocovs, detector = "count", binary.usage = FALSE, trapID = "Detector")

## Check usage
summary(traps)
attr(traps, "usage") # make sure usage looks ok

## Add covariates if desired
covariates(traps)$topography <- traps_df$topography
covariates(traps)$cameratype <- traps_df$cameratype

# (2) Read in the capture history

## Read in
ch_df <- read.csv("data/secr_inputs/ex_capthist.csv") %>% as.data.frame()

# Turn into secr capthist object
capthist <- make.capthist(captures = ch_df, traps = traps) 
summary(capthist, terse = TRUE)
summary(capthist)  

# (3) Make integration mesh/mask
guess_sigma <- 4000
mesh <- make.mask(traps = traps, 
                  buffer = 6 * guess_sigma,
                  spacing = 1500,
                  type = "trapbuffer")

plot(mesh)
plot(traps, add = TRUE)
plot(survey_area, add = T)

## How to choose buffer and spacing? 

## Start with a quick and dirty estimate of sigma
qnd <- autoini(capthist = capthist, mask = make.mask(traps = traps, buffer = 10000, spacing = 1500))
qnd

## Buffer and spacing chosen for "autoini" should not matter much or at all
autoini(capthist = capthist, mask = make.mask(traps = traps, buffer = 10000, spacing = 500))
autoini(capthist = capthist, mask = make.mask(traps = traps, buffer = 50000, spacing = 1500))
autoini(capthist = capthist, mask = make.mask(traps = traps, buffer = 50000, spacing = 5000))

# A reasonable buffer width is 4 * sigma
qnd$sigma * 4

# Suggest.buffer can also suggest a reasonable buffer
suggest.buffer(object = capthist, detectfn = "HN")

# A conservative but reasonable spacing is 0.6 sigma
qnd$sigma * 0.6

## Add mesh covariates from a spatial data source

### Can add from a raster
tri <- raster("data/spatial_covs/TRI.tiff")
tri_crs <- projectRaster(tri, crs = my_crs) # can take a while for big rasters
mesh <- addCovariates(mesh, tri_crs)
names(covariates(mesh))

### Can also add from shapefile
mesh <- addCovariates(mesh, "data/spatial_covs/TRI.shp")
names(covariates(mesh))

### Note TRI in shapefile is at different resolution to raster so values not identical
covariates(mesh)[1:10,1]
covariates(mesh)[1:10,2]

# view model inputs: mask, capthist, traps
plot(mesh) # this is plot.mask
plot(traps(capthist), add = TRUE) # this is plot.traps
plot(capthist, tracks = TRUE, add = TRUE) # this is plot.capthist

# (4) Fit a model

## Start with baseline model, no covariates on anything
slmod.0 <- secr.fit(capthist = capthist, 
                    mask = mesh, 
                    model = list(D~1, g0~1, sigma~1), 
                    detectfn = "HN")

## Examine output
coef(slmod.0) # estimates on link scale
predict(slmod.0) # estimates on natural scale
region.N(slmod.0) # abundance estimates

## Plot detection function
plot(slmod.0, limits = TRUE, xval = seq(0, 30000, length.out = 100))

## Model with density a function of TRI, no covariates on detection
slmod.tri <- secr.fit(capthist = capthist, 
                      mask = mesh, 
                      model = list(D~TRI, g0~1, sigma~1), 
                      detectfn = "HN")

## Examine output
coef(slmod.tri) # estimates on link scale
predict(slmod.tri) # estimates on natural scale
region.N(slmod.tri) # abundance estimates

## Compare models
AIC(slmod.0, slmod.tri)

# (5) Abundance estimate for a different region to the mask used for model fitting
#     (for example, might want to use the survey region)

## Make a new mesh over the region you want to calculate abundance in
newmesh <- make.mask(traps = traps, 
                      buffer = 50000,
                      spacing = 1500,
                      type = "trapbuffer",
                      poly = survey_area)

plot(survey_area)
plot(newmesh, add = TRUE)

## Abundance estimates
region.N(slmod.0, region = newmesh) 

# save(slmod.0, slmod.tri, file = "output/slmods.Rdata")

