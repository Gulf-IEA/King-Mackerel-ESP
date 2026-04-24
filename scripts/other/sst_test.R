
library(abind)
library(dplyr)
library(lubridate)
# library(IEAnalyzeR)
library(here)
library(ggplot2)
library(rerddap)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(ncdf4)
library(cmocean)

### observer data
# setwd("~/data/Fishery_observer_data/ReefExtraction for Coastal Pelagics (i.e. Trollin20240924044830")
# trolling <- read_xlsx('Trolling Reef Observer Program Request (09_24_24).xlsx')
# kmk_tr <- subset(trolling, COMMON_NAME=='MACKEREL, KING' & LAT_BEGIN_SET < 40)
# lons <- range(kmk_tr$LON_BEGIN_SET)
# lats <- range(kmk_tr$LAT_BEGIN_SET)
# obs <- data.frame(lon = kmk_tr$LON_BEGIN_SET,
#                   lat = kmk_tr$LAT_BEGIN_SET)

# define years  --------------------------------
styear <- 1982
enyear <- 2025

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


# download by year to avoid timeout errors --------------------

######################################################
#### don't run while reviewing code; takes awhile ####
#### load saved intermediate files below loop ########
######################################################

review_code <- T ### set to F to rerun download loop

if(review_code == F){
  
  # get ERDDAP info  --------------------------------
  sst <- info('ncdcOisst21Agg_LonPM180') # this may work better
  
  # empty data  -------------------------------------------------
  dat_eez <- dat_swfl <- dat_latx <- dat_desoto <- c()
  
  setwd(here("data/intermediate"))
  system.time(
    for (yr in styear:enyear) { 
      
      cat('\n', yr, '\n')
      
      sst_grab <- griddap(sst, fields = c('anom','sst'), 
                          time = c(paste0(yr,'-01-01'), paste0(yr,'-12-31')),
                          longitude = c(min_lon, max_lon), 
                          latitude = c(min_lat, max_lat), 
                          fmt = 'csv')
      sst_eez_sf <- st_as_sf(sst_grab, coords = c("longitude", "latitude"), crs = 4326)
      
      ### US EEZ
      sst_eez <- st_intersection(sst_eez_sf, gulf_eez2) |>
        st_drop_geometry() |>
        group_by(time) |>
        summarize(sst_degC = mean(sst, na.rm = T),
                  anom_degC = mean(anom, na.rm = T))
      
      ### swfl
      sst_swfl <- st_intersection(sst_eez_sf, swfl) |>
        st_drop_geometry() |>
        # filter(month(time)<3) |> 
        group_by(time) |>
        summarize(sst_degC = mean(sst, na.rm = T),
                  anom_degC = mean(anom, na.rm = T))
      
      ### latx
      sst_latx <- st_intersection(sst_eez_sf, latx) |>
        st_drop_geometry() |>
        # filter(month(time)>5 & month(time)<11) |> 
        group_by(time) |>
        summarize(sst_degC = mean(sst, na.rm = T),
                  anom_degC = mean(anom, na.rm = T))
      
      ### latx
      sst_desoto <- st_intersection(sst_eez_sf, desoto) |>
        st_drop_geometry() |>
        # filter(month(time)>5 & month(time)<11) |> 
        group_by(time) |>
        summarize(sst_degC = mean(sst, na.rm = T),
                  anom_degC = mean(anom, na.rm = T))
      
      
      if (yr == styear) { 
        dat_eez <- sst_eez
        dat_swfl <- sst_swfl
        dat_latx <- sst_latx
        dat_desoto <- sst_desoto
      } 
      else {
        dat_eez <- rbind(dat_eez, sst_eez)
        dat_swfl <- rbind(dat_swfl, sst_swfl)
        dat_latx <- rbind(dat_latx, sst_latx)
        dat_desoto <- rbind(dat_desoto, sst_desoto)
      }
    }
  )
  
  setwd(here("data/intermediate"))
  save(dat_eez, dat_swfl, dat_latx, dat_desoto, file = 'kmk_sst_comb_temp.RData')
  
} else {
  
  setwd(here("data/intermediate"))
  # load('sst_comb_temp2.RData')
  ### load dat_gulf & dat_eez downloaded sst data
  
}


