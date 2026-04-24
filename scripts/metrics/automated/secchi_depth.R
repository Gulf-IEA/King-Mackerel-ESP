
library(abind)
library(lubridate)
library(ncdf4)
library(terra)
library(sf)
library(reticulate)

### area of interest
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
nc_close(bdat)

rbathy <- rast(t(bathy[,ncol(bathy):1]), 
               extent = ext(min_lon, max_lon, min_lat, max_lat), 
               crs = crs)

rbathy[values(rbathy$lyr.1) > (-10)] <- NA
rbathy[values(rbathy$lyr.1) < (-50)] <- NA

mean(values(rbathy$lyr.1), na.rm = T)

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
  st_transform(crs = st_crs(4326)) #|>
# st_simplify(dTolerance = 1000)

wcfl <- crop(gulf_eez2, ext(-86,-81,26,29)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) #|>
# st_simplify(dTolerance = 1000)

desoto <- crop(gulf_eez2, ext(-89,-81,29,31)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) #|>
# st_simplify(dTolerance = 1000)

latx <- crop(gulf_eez2, ext(-99,-89,25,31)) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) #|>
# st_simplify(dTolerance = 1000)

gulf_kmk <- gulf_eez2 |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326)) #|>
# st_simplify(dTolerance = 1000)

gulf_eez <- terra::intersect(eez, iho) |>
  st_as_sf() |> 
  st_transform(crs = st_crs(4326))



# virtualenv_create(envname = "CopernicusMarineR", packages = c("copernicusmarine"))
# Activate the virtual environment (must be done before importing any Python module)
use_virtualenv("CopernicusMarineR", required = TRUE)

# Optional sanity check: confirm which Python reticulate is using
py_config()

# Import the Python module
copernicusmarine <- import("copernicusmarine")

# The adapted command
result <- copernicusmarine$subset(
  dataset_id = "cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M",
  dataset_version="202603",
  variables = list("ZSD"),  # Use list() so reticulate passes a proper Python list
  minimum_longitude = min_lon,
  maximum_longitude = max_lon,
  minimum_latitude = min_lat,
  maximum_latitude = max_lat,
  start_datetime = "1997-09-01T00:00:00",
  end_datetime   = "2025-12-01T00:00:00",
  output_directory = "C:/Users/brendan.turley/Documents/data/copernicusmarine/zsd"
)

setwd('C:/Users/brendan.turley/Documents/data/copernicusmarine/zsd')
dat <- nc_open('cmems_obs-oc_glo_bgc-transp_my_l4-multi-4km_P1M_ZSD_97.98W-80.02W_18.02N-30.98N_1997-09-01-2025-12-01.nc')

zsd <- ncvar_get(dat, 'ZSD')
time <- ncvar_get(dat, 'time') |>
  as.Date(origin = '1900-01-01')
nc_close(dat)

zsd <- aperm(zsd, c(2,1,3))
zsd_r <- rast(zsd[dim(zsd)[1]:1,,], crs="EPSG:4326")
ext(zsd_r) <- c(min_lon, max_lon, min_lat, max_lat)
time(zsd_r) <- as.Date(time)


# a. annual US Gulf EEZ
year_index <- year(time(zsd_r))

ann_gwide <- crop(zsd_r, gulf_kmk) |> mask(gulf_kmk)
ann_gwide_layers <- tapp(ann_gwide, index = year_index, fun = mean, na.rm = TRUE)
ann_gwide_ts <- global(ann_gwide_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

ann_gwide_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  zsd = ann_gwide_ts$mean,
  min = ann_gwide_ts$min,
  max = ann_gwide_ts$max
)

plot(ann_gwide_ts$year, ann_gwide_ts$zsd, typ = 'o', pch = 16,
     panel.first = grid())


# b. winter FL Keys to Naples; also winter off LA
month_index <- month(time(zsd_r))
winter_ind <- which(month_index %in% c(1,2,3))
win_yrs <- year_index[winter_ind]

win_swfl <- crop(zsd_r, swfl) |> mask(swfl)
win_swfl <- win_swfl[[winter_ind]]
win_swfl_layers <- tapp(win_swfl, index=win_yrs, fun=mean, na.rm=TRUE)
win_swfl_ts <- global(win_swfl_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

win_swfl_ts <- data.frame(
  year = unique(year_index)[-1], # Convert index back to Date
  zsd = win_swfl_ts$mean,
  min = win_swfl_ts$min,
  max = win_swfl_ts$max
)

plot(win_swfl_ts$year, win_swfl_ts$zsd, typ = 'o', pch = 16,
     panel.first = grid())


# c. spring (April) and fall (Nov) off Tampa Bay
spring_ind <- which(month_index %in% c(4,5))
spr_yrs <- year_index[spring_ind]

spr_wcfl <- crop(zsd_r, wcfl) |> mask(wcfl)
spr_wcfl <- spr_wcfl[[spring_ind]]
spr_wcfl_layers <- tapp(spr_wcfl, index=spr_yrs, fun=mean, na.rm=TRUE)
spr_wcfl_ts <- global(spr_wcfl_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

spr_wcfl_ts <- data.frame(
  year = unique(year_index)[-1], # Convert index back to Date
  zsd = spr_wcfl_ts$mean,
  min = spr_wcfl_ts$min,
  max = spr_wcfl_ts$max
)

plot(spr_wcfl_ts$year, spr_wcfl_ts$zsd, typ = 'o', pch = 16,
     panel.first = grid())


# d. summer/fall texas to FL panhandle
sufa_ind <- which(month_index %in% c(6,7,8,9,10))
sufa_yrs <- year_index[sufa_ind]

sufa_latx <- crop(zsd_r, latx) |> mask(latx)
sufa_latx <- sufa_latx[[sufa_ind]]
sufa_latx_layers <- tapp(sufa_latx, index=sufa_yrs, fun=mean, na.rm=TRUE)
sufa_latx_ts <- global(sufa_latx_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

sufa_latx_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  zsd = sufa_latx_ts$mean,
  min = sufa_latx_ts$min,
  max = sufa_latx_ts$max
)

plot(sufa_latx_ts$year, sufa_latx_ts$zsd, typ = 'o', pch = 16,
     panel.first = grid())


### plots

mod1 <- lm(zsd~year, data = ann_gwide_ts)
mod1p <- summary(mod1)$coefficients[8]
mod2 <- lm(zsd~year, data = win_swfl_ts)
mod2p <- summary(mod2)$coefficients[8]
mod3 <- lm(zsd~year, data = spr_wcfl_ts)
mod3p <- summary(mod3)$coefficients[8]
mod4 <- lm(zsd~year, data = sufa_latx_ts)
mod4p <- summary(mod4)$coefficients[8]

asp <- NA

setwd("~/R_projects/King-Mackerel-ESP/figures/plots")
png('kmk_zsd.png', width = 8, height = 7, 
    units = 'in', res = 300)
par(mfrow=c(2,2),
    mar = c(3,4,3,1))

plot(ann_gwide_ts$year, ann_gwide_ts$zsd, 
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = 'Secchi Depth (m)',
     main = 'Gulf-wide Annual',
     panel.first = list(if(mod1p<.05){
       abline(mod1, col = 'orange', lwd = 2)
     }))


plot(win_swfl_ts$year, win_swfl_ts$zsd,
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = 'Secchi Depth (m)',
     main = 'SWFL Winter',
     panel.first = list(if(mod2p<.05){
       abline(mod2, col = 'orange', lwd = 2)
     }))

plot(spr_wcfl_ts$year, spr_wcfl_ts$zsd, 
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = 'Secchi Depth (m)',
     main = 'WCFL Spring',
     panel.first = list(if(mod3p<.05){
       abline(mod3, col = 'orange', lwd = 2)
     }))

plot(sufa_latx_ts$year, sufa_latx_ts$zsd,
     typ = 'o', pch = 16, , las = 1, asp = asp,
     xlab = '', ylab = 'Secchi Depth (m)',
     main = 'LATX Summer/Fall',
     panel.first = list(if(mod4p<.05){
       abline(mod4, col = 'orange', lwd = 2)
     }))

dev.off()
