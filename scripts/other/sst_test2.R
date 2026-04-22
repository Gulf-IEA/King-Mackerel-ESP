
library(abind)
library(lubridate)
library(ncdf4)
library(terra)
library(sf)

### this adds a fahrenheit axis on the right of the plot by converting the celcius default
ax_convert_c2f <- function(vals, side = 4, n = 5, las = 1, ...){ ### ... used to pass other parameters for interior fxns
  tick_val <- pretty(vals*(9/5)+32, n = n, ...)
  axis(side, (tick_val-32)*(5/9), tick_val, las = las, ...)
}

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

plot(st_geometry(gulf_kmk))
plot(st_geometry(swfl),add=T,col=2)
plot(st_geometry(wcfl),add=T,col=3)
plot(st_geometry(desoto),add=T,col=4)
plot(st_geometry(latx),add=T,col=5)



setwd("~/R_projects/ESR-indicator-scratch/data/intermediate_files")

for(i in styear:enyear){
  cat(i, '\n')
  tmp <- paste0('sst_',i) |> readRDS()
  
  tmp$sst[which(tmp$sst==-999)] <- NA
  
  if(i==styear){
    sst_a <- tmp$sst
    dates <- tmp$time
  } else {
    sst_a <- abind(sst_a,
                   tmp$sst,
                   along = 3)
    dates <- c(dates,
               tmp$time)
  }
}

# tmp <- paste0('sst_',1994) |> readRDS()
# tmp$sst[which(tmp$sst==-999)] <- NA
# hist(tmp$sst)
# apply(tmp$sst, c(1,2), mean) |> hist()
# apply(tmp$sst, 3, mean,na.rm=T) |> plot()
# tmp$sst[which(tmp$sst>100)]
# which(tmp$sst>100,arr.ind=T)[,3] |> unique()
# hist(tmp$sst[,,98])

sst_a <- aperm(sst_a, c(2,1,3))
sst_r <- rast(sst_a[dim(sst_a)[1]:1,,], crs="EPSG:4326") 
ext(sst_r) <- c(min_lon, max_lon, min_lat, max_lat)
time(sst_r) <- as.Date(dates)

sst_a[which(!is.finite(sst_a))] |> unique()

# for(i in styear:enyear){
#   cat(i, '\n')
#   tmp <- paste0('anom_',i) |> readRDS()
#   
#   if(i==styear){
#     anom_a <- tmp$anom
#     dates <- tmp$time
#   } else {
#     anom_a <- abind(anom_a,
#                     tmp$anom,
#                     along = 3)
#     dates <- c(dates,
#                tmp$time)
#   }
# }
# 
# anom_a <- aperm(anom_a, c(2,1,3))
# anom_r <- rast(anom_a[dim(anom_a)[1]:1,,], crs="EPSG:4326") 
# ext(anom_r) <- c(min_lon, max_lon, min_lat, max_lat)
# time(anom_r) <- as.Date(dates)


# a. annual US Gulf EEZ
year_index <- year(time(sst_r))

ann_gwide <- crop(sst_r, gulf_kmk) |> mask(gulf_kmk)
ann_gwide_layers <- tapp(ann_gwide, index = year_index, fun = mean, na.rm = TRUE)
ann_gwide_ts <- global(ann_gwide_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

ann_gwide_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  ssta = ann_gwide_ts$mean,
  min = ann_gwide_ts$min,
  max = ann_gwide_ts$max
)

plot(ann_gwide_ts$year, ann_gwide_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(ann_gwide_ts$year, rev(ann_gwide_ts$year)), 
                                c(ann_gwide_ts$min, rev(ann_gwide_ts$max)), col = 5),
                        grid()),
     ylim = range(ann_gwide_ts[,2:4]))


ann_gwide_min <- tapp(ann_gwide, index = year_index, fun = min, na.rm = TRUE)
ann_gwide_max <- tapp(ann_gwide, index = year_index, fun = max, na.rm = TRUE)

ann_gwide_min_ts <- global(ann_gwide_min, fun = min, na.rm = TRUE, weighted = T)
ann_gwide_max_ts <- global(ann_gwide_max, fun = max, na.rm = TRUE, weighted = T)

# bday()
gwide_ts <- global(ann_gwide, fun = c('mean',"range"), na.rm = TRUE, weighted = T)
gwide_ts <- data.frame(
  date = time(sst_r), # Convert index back to Date
  ssta = gwide_ts$mean,
  min = gwide_ts$min,
  max = gwide_ts$max
)

aggregate(ssta ~ year(date), data = gwide_ts, 
          function(x) length(which(x>30))) |> plot(typ = 'h')

# b. winter FL Keys to Naples; also winter off LA
month_index <- month(time(sst_r))
winter_ind <- which(month_index %in% c(1,2,3))
win_yrs <- year_index[winter_ind]

win_swfl <- crop(sst_r, swfl) |> mask(swfl)
win_swfl <- win_swfl[[winter_ind]]
win_swfl_layers <- tapp(win_swfl, index=win_yrs, fun=mean, na.rm=TRUE)
win_swfl_ts <- global(win_swfl_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

win_swfl_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  ssta = win_swfl_ts$mean,
  min = win_swfl_ts$min,
  max = win_swfl_ts$max
)

plot(win_swfl_ts$year, win_swfl_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(win_swfl_ts$year, rev(win_swfl_ts$year)), 
                                c(win_swfl_ts$min, rev(win_swfl_ts$max)), col = 4),
                        grid()),
     ylim = range(win_swfl_ts[,2:4]))


win_swfl_min <- tapp(win_swfl, index = win_yrs, fun = min, na.rm = TRUE)
win_swfl_max <- tapp(win_swfl, index = win_yrs, fun = max, na.rm = TRUE)

win_swfl_min_ts <- global(win_swfl_min, fun = min, na.rm = TRUE, weighted = T)
win_swfl_max_ts <- global(win_swfl_max, fun = max, na.rm = TRUE, weighted = T)

plot(win_swfl_ts$year, win_swfl_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(win_swfl_ts$year, rev(win_swfl_ts$year)), 
                                c(win_swfl_min_ts$min, rev(win_swfl_max_ts$max)), col = 4),
                        grid()),
     ylim = range(c(win_swfl_min_ts$min, rev(win_swfl_max_ts$max))))

# c. spring (April) and fall (Nov) off Tampa Bay
spring_ind <- which(month_index %in% c(4,5))
spr_yrs <- year_index[spring_ind]

spr_wcfl <- crop(sst_r, wcfl) |> mask(wcfl)
spr_wcfl <- spr_wcfl[[spring_ind]]
spr_wcfl_layers <- tapp(spr_wcfl, index=spr_yrs, fun=mean, na.rm=TRUE)
spr_wcfl_ts <- global(spr_wcfl_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

spr_wcfl_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  ssta = spr_wcfl_ts$mean,
  min = spr_wcfl_ts$min,
  max = spr_wcfl_ts$max
)

plot(spr_wcfl_ts$year, spr_wcfl_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(spr_wcfl_ts$year, rev(spr_wcfl_ts$year)), 
                                c(spr_wcfl_ts$min, rev(spr_wcfl_ts$max)), col = 3),
                        grid()),
     ylim = range(spr_wcfl_ts[,2:4]))

spr_wcfl_min <- tapp(spr_wcfl, index = spr_yrs, fun = min, na.rm = TRUE)
spr_wcfl_max <- tapp(spr_wcfl, index = spr_yrs, fun = max, na.rm = TRUE)

spr_wcfl_min_ts <- global(spr_wcfl_min, fun = min, na.rm = TRUE, weighted = T)
spr_wcfl_max_ts <- global(spr_wcfl_max, fun = max, na.rm = TRUE, weighted = T)

plot(spr_wcfl_ts$year, spr_wcfl_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(spr_wcfl_ts$year, rev(spr_wcfl_ts$year)), 
                                c(spr_wcfl_min_ts$min, rev(spr_wcfl_max_ts$max)), col = 3),
                        grid()),
     ylim = range(c(spr_wcfl_min_ts$min, rev(spr_wcfl_max_ts$max))))

# d. summer/fall texas to FL panhandle
sufa_ind <- which(month_index %in% c(6,7,8,9,10))
sufa_yrs <- year_index[sufa_ind]

sufa_latx <- crop(sst_r, latx) |> mask(latx)
sufa_latx <- sufa_latx[[sufa_ind]]
sufa_latx_layers <- tapp(sufa_latx, index=sufa_yrs, fun=mean, na.rm=TRUE)
sufa_latx_ts <- global(sufa_latx_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

sufa_latx_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  ssta = sufa_latx_ts$mean,
  min = sufa_latx_ts$min,
  max = sufa_latx_ts$max
)

plot(sufa_latx_ts$year, sufa_latx_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(sufa_latx_ts$year, rev(sufa_latx_ts$year)), 
                                c(sufa_latx_ts$min, rev(sufa_latx_ts$max)), col = 2),
                        grid()),
     ylim = range(sufa_latx_ts[,2:4]))

sufa_latx_min <- tapp(sufa_latx, index = sufa_yrs, fun = min, na.rm = TRUE)
sufa_latx_max <- tapp(sufa_latx, index = sufa_yrs, fun = max, na.rm = TRUE)

sufa_latx_min_ts <- global(sufa_latx_min, fun = min, na.rm = TRUE, weighted = T)
sufa_latx_max_ts <- global(sufa_latx_max, fun = max, na.rm = TRUE, weighted = T)

plot(sufa_latx_ts$year, sufa_latx_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(sufa_latx_ts$year, rev(sufa_latx_ts$year)), 
                                c(sufa_latx_min_ts$min, rev(sufa_latx_max_ts$max)), col = 3),
                        grid()),
     ylim = range(c(sufa_latx_min_ts$min, rev(sufa_latx_max_ts$max))))


sufa_ds <- crop(sst_r, desoto) |> mask(desoto)
sufa_ds <- sufa_ds[[sufa_ind]]
sufa_ds_layers <- tapp(sufa_ds, index=sufa_yrs, fun=mean, na.rm=TRUE)
sufa_ds_ts <- global(sufa_ds_layers, fun = c('mean',"range"), na.rm = TRUE, weighted = T)

sufa_ds_ts <- data.frame(
  year = unique(year_index), # Convert index back to Date
  ssta = sufa_ds_ts$mean,
  min = sufa_ds_ts$min,
  max = sufa_ds_ts$max
)

plot(sufa_ds_ts$year, sufa_ds_ts$ssta, typ = 'o', pch = 16, 
     panel.first = list(polygon(c(sufa_ds_ts$year, rev(sufa_ds_ts$year)), 
                                c(sufa_ds_ts$min, rev(sufa_ds_ts$max)), col = 2),
                        grid()),
     ylim = range(sufa_ds_ts[,2:4]))



# sufa_latx_ts <- global(sufa_latx_layers, fun = quantile, probs = c(0,.5,1), na.rm = TRUE, weighted = T)
# 
# sufa_latx_ts <- data.frame(
#   year = unique(year_index), # Convert index back to Date
#   ssta = sufa_latx_ts$X50.,
#   min = sufa_latx_ts$X0.,
#   max = sufa_latx_ts$X100.
# )
# 
# plot(sufa_latx_ts$year, sufa_latx_ts$ssta, typ = 'o', pch = 16, 
#      panel.first = list(polygon(c(sufa_latx_ts$year, rev(sufa_latx_ts$year)), 
#                                 c(sufa_latx_ts$min, rev(sufa_latx_ts$max)), col = 2),
#                         grid()),
#      ylim = range(sufa_latx_ts[,2:4]))


mod1 <- lm(ssta~year, data = ann_gwide_ts)
mod1p <- summary(mod1)$coefficients[8]
mod2 <- lm(ssta~year, data = win_swfl_ts)
mod2p <- summary(mod2)$coefficients[8]
mod3 <- lm(ssta~year, data = spr_wcfl_ts)
mod3p <- summary(mod3)$coefficients[8]
mod4 <- lm(ssta~year, data = sufa_latx_ts)
mod4p <- summary(mod4)$coefficients[8]

asp <- 9/1

setwd("~/R_projects/King-Mackerel-ESP/figures/plots")
png('kmk_sst.png', width = 8, height = 7, 
    units = 'in', res = 300)
par(mfrow=c(2,2),
    mar = c(3,3,3,3))

plot(ann_gwide_ts$year, ann_gwide_ts$ssta, 
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = '',
     main = 'Gulf-wide Annual',
     panel.first = list(if(mod1p<.05){
       abline(mod1, col = 'orange', lwd = 2)
     }))
mtext('(°C)', side = 3, adj = -.1, line = .5)
mtext('(°F)', side = 3, adj = 1.1, line = .5)
ax_convert_c2f(axTicks(side = 2), n = 4)


plot(win_swfl_ts$year, win_swfl_ts$ssta,
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = '',
     main = 'SWFL Winter',
     panel.first = list(if(mod2p<.05){
       abline(mod2, col = 'orange', lwd = 2)
     }))
mtext('(°C)', side = 3, adj = -.1, line = .5)
mtext('(°F)', side = 3, adj = 1.1, line = .5)
ax_convert_c2f(axTicks(side = 2), n = 4)

plot(spr_wcfl_ts$year, spr_wcfl_ts$ssta, 
     typ = 'o', pch = 16, las = 1, asp = asp,
     xlab = '', ylab = '',
     main = 'WCFL Spring',
     panel.first = list(if(mod3p<.05){
       abline(mod3, col = 'orange', lwd = 2)
     }))
mtext('(°C)', side = 3, adj = -.1, line = .5)
mtext('(°F)', side = 3, adj = 1.1, line = .5)
ax_convert_c2f(axTicks(side = 2), n = 4)

plot(sufa_latx_ts$year, sufa_latx_ts$ssta,
     typ = 'o', pch = 16, , las = 1, asp = asp,
     xlab = '', ylab = '',
     main = 'LATX Summer/Fall',
     panel.first = list(if(mod4p<.05){
       abline(mod4, col = 'orange', lwd = 2)
     }))
mtext('(°C)', side = 3, adj = -.1, line = .5)
mtext('(°F)', side = 3, adj = 1.1, line = .5)
ax_convert_c2f(axTicks(side = 2), n = 4)

dev.off()




