
library(ncdf4)
library(terra)
library(sf)


# define spatial domain  --------------------------------
min_lon <- -98
max_lon <- -80
min_lat <- 18
max_lat <- 31

### bathymetry

burl <- 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc'

bdat <- nc_open(burl)
crs <- 'EPSG:4326'

ln <- ncvar_get(bdat, 'lon')
ln_i <- which(ln>=min_lon & ln<=max_lon)
lt <- ncvar_get(bdat, 'lat')
lt_i <- which(lt>=min_lat & lt<=max_lat)

bathy <- ncvar_get(bdat, 'z', 
                   start = c(ln_i[1],lt_i[1]),
                   count = c(length(ln_i), length(lt_i)))

# blon <- ln[ln_i]
# blat <- lt[lt_i]

rbathy <- rast(t(bathy[,ncol(bathy):1]), 
               extent = ext(min_lon, max_lon, min_lat, max_lat), 
               crs = crs)
# rbathy2 <- rbathy
# rbathy2[values(rbathy2$lyr.1) > 0] <- NA
# rbathy2[values(rbathy2$lyr.1) < (-200)] <- NA
# plot(rbathy2)

# depths <- terra::extract(rbathy, obs)
# depths$lyr.1[depths$lyr.1>0] <- NA
# hist(depths$lyr.1)
# summary(depths$lyr.1)
# quantile(depths$lyr.1, c(.025,.975), na.rm = T)
# ### isolate shelf
# rbathy[values(rbathy$lyr.1) > (-10)] <- NA
# rbathy[values(rbathy$lyr.1) < (-40)] <- NA

rbathy[values(rbathy$lyr.1) > (-10)] <- NA
rbathy[values(rbathy$lyr.1) < (-50)] <- NA

plot(rbathy,
     main = 'Continential Shelf cutout')

### set all values to 1; then make into shapefile
rbathy[!is.na(values(rbathy$lyr.1))] <- 1
rbathy_p <- as.polygons(rbathy) |> makeValid()

# load shapefile to subset  --------------------------------
### shapefiles downloaded from marineregions.org (future goal implement mregions2 R package for shapefile)
setwd("~/data/shapefiles/gulf_eez")
eez <- vect('eez.shp') |> makeValid()

setwd("~/data/shapefiles/gulf_iho")
iho <- vect('iho.shp') |> makeValid()

gulf_eez <- terra::intersect(eez, iho)

gulf_eez2 <- terra::intersect(gulf_eez, rbathy_p)
swfl <- crop(gulf_eez2, ext(-84,-81,24,26)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) |>
  st_simplify(dTolerance = 1000)
desoto <- crop(gulf_eez2, ext(-89,-85,27,31)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) |>
  st_simplify(dTolerance = 1000)
latx <- crop(gulf_eez2, ext(-96,-89,27,31)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) |>
  st_simplify(dTolerance = 1000)

gulf_eez2 <- gulf_eez2 |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) |>
  st_simplify(dTolerance = 1000)



setwd("~/R_projects/ESR-indicator-scratch/data/intermediate_files")

anom_1982 <- readRDS('anom_1982')
anom_1982_a <- aperm(anom_1982$anom, c(2,1,3))
anom82_r <- rast(anom_1982_a[dim(anom_1982_a)[1]:1,,], crs="EPSG:4326") 
ext(anom82_r) <- c(min_lon, max_lon, min_lat, max_lat)
time(anom82_r) <- as.Date(anom_1982$time)
month_index <- format(time(anom82_r), "%Y-%m")

sst_subset <- crop(anom82_r, gulf_eez2) |> mask(gulf_eez2)
sst_monthly_layers <- tapp(sst_subset, index = month_index, fun = mean, na.rm = TRUE)
ts_values <- global(sst_monthly_layers, fun = "mean", na.rm = TRUE, weighted = T)

sst_timeseries <- data.frame(
  Date = as.Date(paste0(unique(month_index), "-01")), # Convert index back to Date
  SST = ts_values$mean
)


anom_1982 <- readRDS('anom_1982')
anom_1983 <- readRDS('anom_1983')

test <- abind(anom_1982$anom,
                anom_1983$anom,
                along=3)
dates <- c(anom_1982$time,
           anom_1983$time)

test2 <- aperm(test, c(2,1,3))
test_r <- rast(test2[dim(test2)[1]:1,,], crs="EPSG:4326") 
ext(test_r) <- c(min_lon, max_lon, min_lat, max_lat)
time(test_r) <- as.Date(dates)
month_index <- format(time(test_r), "%Y-%m")

sst_subset <- crop(test_r, gulf_eez2) |> mask(gulf_eez2)
sst_monthly_layers <- tapp(sst_subset, index = month_index, fun = mean, na.rm = TRUE)
ts_values <- global(sst_monthly_layers, fun = "mean", na.rm = TRUE, weighted = T)

sst_timeseries <- data.frame(
  Date = as.Date(paste0(unique(month_index), "-01")), # Convert index back to Date
  SST = ts_values$mean
)


### one way to check this would be to recreate the Gulf ESR monthly index and compare

