
library(dplyr)
library(here)
library(readxl)
library(sf)
library(scatterpie)
library(ggplot2)
library(terra)

setwd("C:/Users/brendan.turley/Documents/data/shapefiles/cflp_statgrid")
sz_shp <- vect('CFLP_StatGrid_2013_v20140210.shp') |>
  st_as_sf()
sz_shp$AREA_FISHED <- sz_shp$SZ_ID

setwd("~/data/shapefiles/ne_10m_admin_1_states_provinces")
states <- vect('ne_10m_admin_1_states_provinces.shp') |> makeValid() |>
  st_as_sf() 
states <- st_crop(states, xmin = -100, ymin = 24, xmax = -80, ymax = 31)


### boem platforms
setwd("~/R_projects/ESR-indicator-scratch/data/processed")
plat <- read.csv('platforms_yr.csv')

png(here(paste0("figures/plots/platfroms_plot.png")),
    width = 10, height = 4, units = 'in', res = 300)
b1 <- barplot(plat$nplt, 
              col = 'gray50', las = 2,
              xlab = '', ylab = 'Number of platforms')
axis(1, b1[seq(4,length(b1),5)], seq(1945, 2025, 5))
dev.off()


### hypoxia extent
setwd("~/R_projects/Gulf-ESR/data/unformatted")
hyp <- read_xlsx('Historical_Hypoxia_Extent_data_2025.xlsx')
hyp$`NOAA Figure` <- as.numeric(hyp$`NOAA Figure`)
hyp <- hyp[, c('Year','NOAA Figure')] |> 
  setNames(c('year','hypoxia_extent_mi2')) |>
  type.convert()

png(here(paste0("figures/plots/hypoxia_plot.png")),
    width = 7, height = 4, units = 'in', res = 300)
b2 <- barplot(hyp$hypoxia_extent_mi2, 
              col = 'orangered3', las = 2,
              xlab = '', ylab = 'Hypoxia Extent (square miles)')
axis(1, b2[seq(1,length(b2),5)], seq(1985, 2025, 5))
dev.off()


### cflp

#### read data and subset ####--------------------------------------------------
gom_st <- c('FL', 'AL', 'MS', 'LA', 'TX')

setwd("C:/Users/brendan.turley/Documents/CMP/data/cflp")
cflp <- readRDS('CFLPblake.rds')
# cflp_ne <- subset(cflp, LAND_YEAR<2024 & LAND_YEAR>1999 & CATCH_TYPE == 'CATCH') |>
#   subset(REGION == 'NATL' & is.element(ST_ABRV, natl_st))
cflp <- subset(cflp, LAND_YEAR>1999 & CATCH_TYPE == 'CATCH') |>
  subset(REGION == 'GOM' & is.element(ST_ABRV, gom_st))
gc()

### pull out handlines only
# table(cflp$GEAR)
gear_keep <- c('H', 'E', 'TR')
# gear_keep <- c('TR')
cflp_hl <- subset(cflp , is.element(cflp$GEAR, gear_keep)) |>
  subset(FLAG_MULTIGEAR==0 & FLAG_MULTIAREA==0)
# saveRDS(cflp_hl, 'cflp_gulfsa_temp.rds')
## cflp_hl <- readRDS('cflp_gulfsa_temp.rds')
# gc()

#### CPUE calculation and days away correction ####-----------------------------
## following methods by Walter & McCarthy 2014 (1993-2013SEDAR38-DW-10)
## CPUE = total kilograms of king mackerel/(number of lines fished*number of hooks per line*total hours fished)
## total whole pounds seems most appropriate; other 2 have lots of zeros
cflp_hl$tot_kg <- cflp_hl$TOTAL_WHOLE_POUNDS / 2.205
cflp_hl$pue <- (cflp_hl$NUMGEAR * cflp_hl$EFFORT * cflp_hl$FISHED)
cflp_hl$cpue <- cflp_hl$tot_kg / cflp_hl$pue # kg catch per number of hooks X hours
cflp_hl$cpue <- ifelse(is.infinite(cflp_hl$cpue), NA, cflp_hl$cpue)

### correct days_away
diff_time <- (cflp_hl$LAND_DATE - cflp_hl$DEPART_DATE)
units(diff_time) <- 'days'
cflp_hl$days_away_corrected <- round(diff_time + 1)
rm(diff_time)
gc()
### end


#### pull vessels landing the top 80% of KMK ####-------------------------------
### this sorts by region and state
### Walter sorts by year and then total landings to get at the 80%
vessel_yrs <- with(subset(cflp_hl, COMMON_NAME=='MACKERELS, KING AND CERO',
                          select = c(LAND_YEAR, VESSEL_ID, REGION, ST_ABRV)),
                   aggregate(LAND_YEAR ~ VESSEL_ID + REGION + ST_ABRV,
                             FUN = function(x) length(unique(x))))
vessel_tot <- with(subset(cflp_hl,
                          select = c(tot_kg, VESSEL_ID, REGION, ST_ABRV)),
                   aggregate(tot_kg ~ VESSEL_ID + REGION + ST_ABRV, FUN = sum, na.rm = T))
vessel_kmk_tot <- with(subset(cflp_hl, COMMON_NAME=='MACKERELS, KING AND CERO',
                              select = c(tot_kg, VESSEL_ID, REGION, ST_ABRV)),
                       aggregate(tot_kg ~ VESSEL_ID + REGION + ST_ABRV, FUN = sum, na.rm = T)) |>
  setNames(c("VESSEL_ID","REGION",'ST_ABRV',"kmk_tot_kg"))
vessel_select <- merge(vessel_yrs, vessel_tot,
                       by = c('VESSEL_ID', 'REGION','ST_ABRV')) |>
  merge(vessel_kmk_tot, by = c('VESSEL_ID','REGION','ST_ABRV'))
vessel_select$kmk_pro <- vessel_select$kmk_tot_kg / vessel_select$tot_kg
vessel_select <- vessel_select[order(vessel_select$LAND_YEAR,
                                     vessel_select$kmk_tot_kg,
                                     vessel_select$kmk_pro,
                                     decreasing = T), ]
n <- 1
ves_id <- list()
for(i in c('GOM')){
  ti <- subset(vessel_select, REGION==i)
  if(i=='GOM') st_sl <- gom_st
  # if(i=='SATL') st_sl <- satl_st
  for(j in st_sl){
    tj <- subset(ti, ST_ABRV==j)
    tj$cummulative <- cumsum(tj$tot_kg)/sum(tj$tot_kg)
    ves_id[[n]] <- tj$VESSEL_ID[which(tj$cummulative<=.8)]
    n <- n + 1
  }
}
kmk_ves <- unique(unlist(ves_id))

cflp_hl_0 <- cflp_hl[is.element(cflp_hl$VESSEL_ID, kmk_ves), ] |>
  subset(
    NUMGEAR < quantile(cflp_hl$NUMGEAR, .995, na.rm = T) &
      EFFORT < quantile(cflp_hl$EFFORT, .995, na.rm = T) &
      FISHED < quantile(cflp_hl$FISHED, .995, na.rm = T) &
      tot_kg < quantile(cflp_hl$tot_kg, .995, na.rm = T) &
      days_away_corrected < quantile(cflp_hl$days_away_corrected, .995, na.rm = T)
  )
cflp_hl_1 <- cflp_hl_0
rm(cflp_hl, cflp_hl_0)
gc()

#### cpue per region overtime ####----------------------------------------------

## CPUE = total kilograms of king mackerel/(number of lines fished*number of hooks per line*total hours fished)
cpue_yr <- aggregate(cbind(cpue, pue, NUMGEAR, EFFORT, FISHED, tot_kg) ~ LAND_YEAR,
                     data = subset(cflp_hl_1, 
                                   COMMON_NAME=='MACKERELS, KING AND CERO' &
                                     REGION=='GOM'),
                     mean, na.rm = T)
tot_yr <- aggregate(cbind(FISHED, tot_kg) ~ LAND_YEAR,
                    data = subset(cflp_hl_1, 
                                  COMMON_NAME=='MACKERELS, KING AND CERO' &
                                    REGION=='GOM'),
                    sum, na.rm = T)

labels <- c('Year', 
            'Mean CPUE (kg / hook-hours)',
            'Mean Effort (hook hours)',
            "Mean Number of lines fished",
            'Mean Number of hooks per line',
            'Mean hours fished',
            'Mean landings (kg)')

for(i in 2:7){
  png(here(paste0("figures/plots/mean_",names(cpue_yr)[i],"_plot.png")),
      width = 7, height = 4, units = 'in', res = 300)
  plot(cpue_yr$LAND_YEAR, cpue_yr[,i],
       typ = 'o', pch = 16, las = 1,
       xlab = 'Year', ylab = labels[i])
  abline(h = mean(cpue_yr[,i]),
         lty = 2)
  grid()
  abline(lm(cpue_yr[,i] ~ cpue_yr$LAND_YEAR), col = 1, lwd = 2)
  dev.off()
}

labels2 <- c('Year',
            'Total hours fished',
            'Total landings (kg)')

for(i in 2:3){
  png(here(paste0("figures/plots/tot_",names(tot_yr)[i],"_plot.png")),
      width = 7, height = 4, units = 'in', res = 300)
  plot(tot_yr$LAND_YEAR, tot_yr[,i],
       typ = 'o', pch = 16,
       xlab = 'Year', ylab = labels2[i])
  abline(h = mean(tot_yr[,i]),
         lty = 2)
  grid()
  abline(lm(tot_yr[,i] ~ tot_yr$LAND_YEAR), col = 1, lwd = 2)
  dev.off()
}

### end



#### cpue per state overtime ####----------------------------------------------

## CPUE = total kilograms of king mackerel/(number of lines fished*number of hooks per line*total hours fished)
cpue_st_yr <- aggregate(cbind(cpue, pue, NUMGEAR, EFFORT, FISHED, tot_kg) ~ LAND_YEAR + ST_ABRV,
                        data = subset(cflp_hl_1, 
                                      COMMON_NAME=='MACKERELS, KING AND CERO' &
                                        REGION=='GOM'),
                        mean, na.rm = T)
tot_st_yr <- aggregate(cbind(FISHED, tot_kg) ~ LAND_YEAR + ST_ABRV,
                       data = subset(cflp_hl_1, 
                                     COMMON_NAME=='MACKERELS, KING AND CERO' &
                                       REGION=='GOM'),
                       sum, na.rm = T)

labels <- c('Year', 
            'State',
            'Mean CPUE (kg / hook-hours)',
            'Mean Effort (hook hours)',
            "Mean Number of lines fished",
            'Mean Number of hooks per line',
            'Mean hours fished',
            'Mean landings (kg)')

gulf <- c('AL', 'FL', 'LA', 'MS', 'TX')
for(i in 3:8){
  png(here(paste0("figures/plots/mean_st_",names(cpue_st_yr)[i],"_plot.png")),
      width = 7, height = 4, units = 'in', res = 300)
  plot(cpue_st_yr$LAND_YEAR, cpue_st_yr[,i], typ = 'n', xlab = 'Year', ylab = labels[i])
  for(j in 1:length(gulf)){
    tmp <- subset(cpue_st_yr, ST_ABRV==gulf[j])
    points(tmp$LAND_YEAR, tmp[,i], typ = 'o', col = j, pch = 16, lwd = 2)
  }
  grid()
  legend('topleft',gulf, lty=1, col=1:length(gulf),bty = 'n', pch = 16, lwd = 2)
  dev.off()
}

labels2 <- c('Year',
             'State',
             'Total hours fished',
             'Total landings (kg)')

for(i in 3:4){
  png(here(paste0("figures/plots/tot_st_",names(tot_st_yr)[i],"_plot.png")),
      width = 7, height = 4, units = 'in', res = 300)
  plot(tot_st_yr$LAND_YEAR, tot_st_yr[,i], typ = 'n', xlab = 'year', ylab = labels2[i])
  for(j in 1:length(gulf)){
    tmp <- subset(tot_st_yr, ST_ABRV==gulf[j])
    points(tmp$LAND_YEAR, tmp[,i], typ = 'o', col = j, pch = 16, lwd = 2)
  }
  grid()
  legend('topleft',gulf, lty=1, col=1:length(gulf),bty = 'n', pch = 16, lwd = 2)
  dev.off()
}

### end



### has the length of trips changed or have there been less days fished?
days_yr <- aggregate(days_away_corrected ~ LAND_YEAR,
                     data = subset(cflp_hl_1, 
                                   COMMON_NAME=='MACKERELS, KING AND CERO' &
                                     REGION=='GOM'),
                     mean, na.rm = T) #median value not informative

days_yr_st <- aggregate(days_away_corrected ~ LAND_YEAR + ST_ABRV,
                        data = subset(cflp_hl_1, 
                                      COMMON_NAME=='MACKERELS, KING AND CERO' &
                                        REGION=='GOM'),
                        mean, na.rm = T) #median value not informative

trps_yr <- aggregate(SCHEDULE_NUMBER ~ LAND_YEAR,
                     data = subset(cflp_hl_1, 
                                   COMMON_NAME=='MACKERELS, KING AND CERO' &
                                     REGION=='GOM'),
                     function(x) length(unique(x)))

trps_yr_st <- aggregate(SCHEDULE_NUMBER ~ LAND_YEAR + ST_ABRV,
                     data = subset(cflp_hl_1, 
                                   COMMON_NAME=='MACKERELS, KING AND CERO' &
                                     REGION=='GOM'),
                     function(x) length(unique(x)))


png(here(paste0("figures/plots/mean_lth_trip_plot.png")),
    width = 7, height = 4, units = 'in', res = 300)
plot(days_yr$LAND_YEAR, days_yr$days_away_corrected,
     typ = 'o', pch = 16,
     xlab = 'year', ylab = 'Length of trip (days)')
dev.off()

png(here(paste0("figures/plots/mean_st_lth_trip_plot.png")),
    width = 7, height = 4, units = 'in', res = 300)
plot(days_yr_st$LAND_YEAR, days_yr_st$days_away_corrected,
     typ = 'n',
     xlab = 'year', ylab = 'Length of trip (days)')
gulf <- c('AL', 'FL', 'LA', 'MS', 'TX')
# gulf <- c('AL', 'FL', 'LA', 'MS')
for(j in 1:length(gulf)){
  tmp <- subset(days_yr_st, ST_ABRV==gulf[j])
  points(tmp$LAND_YEAR, tmp$days_away_corrected, typ = 'o', col = j, pch = 16, lwd = 2)
}
legend('topleft',gulf, lty=1, col=1:length(gulf),bty = 'n', pch = 16, lwd = 2)
grid()
dev.off()


png(here(paste0("figures/plots/mean_num_trip_plot.png")),
    width = 7, height = 4, units = 'in', res = 300)
plot(trps_yr$LAND_YEAR, trps_yr$SCHEDULE_NUMBER,
     typ = 'o', pch = 16,
     xlab = 'year', ylab = 'Number of Trips')
dev.off()


png(here(paste0("figures/plots/mean_st_num_trip_plot.png")),
    width = 7, height = 4, units = 'in', res = 300)
plot(trps_yr_st$LAND_YEAR, trps_yr_st$SCHEDULE_NUMBER,
     typ = 'n',
     xlab = 'year', ylab = 'Number of Trips')
gulf <- c('AL', 'FL', 'LA', 'MS', 'TX')
# gulf <- c('AL', 'FL', 'LA', 'MS')
for(j in 1:length(gulf)){
  tmp <- subset(trps_yr_st, ST_ABRV==gulf[j])
  points(tmp$LAND_YEAR, tmp$SCHEDULE_NUMBER, typ = 'o', col = j, pch = 16, lwd = 2)
}
legend('topleft',gulf, lty=1, col=1:length(gulf),bty = 'n', pch = 16, lwd = 2)
grid()
dev.off()



#### hovmoller plots month x year proportion of total catch per month per year ####----------------------------------

# aggregrate by state, month, year

yr_mth <- expand.grid(LAND_MONTH = 1:12, LAND_YEAR = 2000:2024)

cpue_hov1 <- subset(cflp_hl_1, COMMON_NAME=='MACKERELS, KING AND CERO') |>
  aggregate(cpue ~ LAND_YEAR + LAND_MONTH + ST_ABRV + REGION,
            mean, na.rm = T)

gulf_fl_cpue <- subset(cpue_hov1, ST_ABRV=='FL' & REGION=='GOM') |>
  merge(yr_mth, by = c('LAND_YEAR', 'LAND_MONTH'), all = T) |>
  arrange(LAND_YEAR, LAND_MONTH)

gulf_la_cpue <- subset(cpue_hov1, ST_ABRV=='LA' & REGION=='GOM') |>
  merge(yr_mth, by = c('LAND_YEAR', 'LAND_MONTH'), all = T) |>
  arrange(LAND_YEAR, LAND_MONTH)

image(matrix(gulf_fl_cpue$cpue, 25, 12, byrow = T))
image(matrix(gulf_la_cpue$cpue, 25, 12, byrow = T))


landings_hov1 <- subset(cflp_hl_1, COMMON_NAME=='MACKERELS, KING AND CERO') |>
  aggregate(tot_kg ~ LAND_YEAR + LAND_MONTH + ST_ABRV + REGION,
            sum, na.rm = T)

gulf_fl_landings <- subset(landings_hov1, ST_ABRV=='FL' & REGION=='GOM') |>
  merge(yr_mth, by = c('LAND_YEAR', 'LAND_MONTH'), all = T) |>
  arrange(LAND_YEAR, LAND_MONTH)
gfl_hov_landings <- matrix(gulf_fl_landings$tot_kg, 25, 12, byrow = T)
gfl_hov_plandings <- t(t(gfl_hov_landings) / apply(gfl_hov_landings, 1, sum, na.rm = T))

gulf_la_landings <- subset(landings_hov1, ST_ABRV=='LA' & REGION=='GOM') |>
  merge(yr_mth, by = c('LAND_YEAR', 'LAND_MONTH'), all = T) |>
  arrange(LAND_YEAR, LAND_MONTH)
gla_hov_landings <- matrix(gulf_la_landings$tot_kg, 25, 12, byrow = T)
gla_hov_plandings <- t(t(gla_hov_landings) / apply(gla_hov_landings, 1, sum, na.rm = T))

gulf_tx_landings <- subset(landings_hov1, ST_ABRV=='TX' & REGION=='GOM') |>
  merge(yr_mth, by = c('LAND_YEAR', 'LAND_MONTH'), all = T) |>
  arrange(LAND_YEAR, LAND_MONTH)
gtx_hov_landings <- matrix(gulf_tx_landings$tot_kg, 25, 12, byrow = T)
gtx_hov_plandings <- t(t(gtx_hov_landings) / apply(gtx_hov_landings, 1, sum, na.rm = T))

image(2000:2024, 1:12,
      gfl_hov_plandings,
      las = 1, xlab = 'year', ylab = 'month',
      breaks = seq(0,.8,.05), col = cmocean('dense')(16))

image(2000:2024, 1:12,
      gla_hov_plandings,
      las = 1, xlab = 'year', ylab = 'month',
      breaks = seq(0,.8,.05), col = cmocean('dense')(16))

image(2000:2024, 1:12,
      gtx_hov_plandings,
      las = 1, xlab = 'year', ylab = 'month',
      breaks = seq(0,.8,.05), col = cmocean('dense')(16))


b <- boxplot(cpue ~ LAND_YEAR, data = subset(cflp_hl_1,
                                                      COMMON_NAME=='MACKERELS, KING AND CERO'),
             pch = 16, lty = 1, varwidth = F, staplewex = 0, lwd = 2, outline = F)
b <- boxplot(tot_kg ~ LAND_YEAR, data = subset(cflp_hl_1,
                                             COMMON_NAME=='MACKERELS, KING AND CERO'),
             pch = 16, lty = 1, varwidth = F, staplewex = 0, lwd = 2, outline = F)


### where kmk are caught per season
gulf_kmk_trips <- subset(cflp_hl_1, 
                         COMMON_NAME=='MACKERELS, KING AND CERO' & REGION=='GOM')
tot_land <- aggregate(tot_kg ~ SCHEDULE_NUMBER, #+
                      # REGION + ST_ABRV + AREA_FISHED +
                      # LAND_YEAR + LAND_MONTH,
                      data = subset(cflp_hl_1, 
                                    SCHEDULE_NUMBER %in% gulf_kmk_trips$SCHEDULE_NUMBER),
                      sum, na.rm = T)
tot_kmk_land <- aggregate(tot_kg ~ SCHEDULE_NUMBER, #+
                          # REGION + ST_ABRV + AREA_FISHED +
                          # LAND_YEAR + LAND_MONTH,
                          data = subset(cflp_hl_1, 
                                        SCHEDULE_NUMBER %in% gulf_kmk_trips$SCHEDULE_NUMBER &
                                          COMMON_NAME=='MACKERELS, KING AND CERO'),
                          sum, na.rm = T)
names(tot_kmk_land)[length(names(tot_kmk_land))] <- 'tot_kmk_kg'
tot_landm <- merge(tot_land, tot_kmk_land)
tot_landm$kmk_pro <- tot_landm$tot_kmk_kg / tot_landm$tot_kg

kmk_trips <- subset(tot_landm, kmk_pro > .9, select = 'SCHEDULE_NUMBER')
kmk_trips_catch <- subset(cflp_hl_1, 
                          SCHEDULE_NUMBER %in% kmk_trips$SCHEDULE_NUMBER)

kmk_trips_catch <- kmk_trips_catch |>
  mutate(
    season = case_when(
      LAND_MONTH < 3 ~ 'win',
      LAND_MONTH > 2 & LAND_MONTH < 6 ~ 'spr',
      LAND_MONTH > 5 & LAND_MONTH < 9 ~ 'sum',
      LAND_MONTH > 8 & LAND_MONTH < 12 ~ 'aut',
      LAND_MONTH == 12 ~ 'win'))

tot_sea_area <- aggregate(tot_kg ~ season + AREA_FISHED,
                          data = subset(kmk_trips_catch,
                                        REGION=='GOM' &
                                          COMMON_NAME=='MACKERELS, KING AND CERO' &
                                          LAND_YEAR>2012),
                          # mean, na.rm = T)
                          sum, na.rm = T)

dat_sea <- reshape(tot_sea_area,timevar = 'season', idvar = 'AREA_FISHED', direction = 'wide')
dat_sea[is.na(dat_sea)] <- 0

dat_sea_sf <-  merge(dat_sea,
                     sz_shp,
                     by = c('AREA_FISHED')) |>
  st_as_sf()
# dat_sea_sf$centroids <- st_coordinates(st_centroid(dat_sea_sf))
dat_sea_sf$lon <- st_coordinates(st_centroid(dat_sea_sf))[,'X']
dat_sea_sf$lat <- st_coordinates(st_centroid(dat_sea_sf))[,'Y']

test <- dat_sea_sf |> 
  select(lon, lat, tot_kg.win, tot_kg.spr, tot_kg.sum, tot_kg.aut) |>
  st_drop_geometry()

theme_set(theme_bw())

ggplot(data = dat_sea_sf) + geom_sf()  + 
  geom_sf(data = states, fill = 'gray40', color = 1) + 
  coord_sf(xlim = c(-97, -81), ylim = c(24, 31)) +
  geom_scatterpie(aes(x=lon, y=lat),
                  data=test,
                  cols=c('tot_kg.win', 'tot_kg.spr', 'tot_kg.sum', 'tot_kg.aut'),
                  color = 1) +
  scale_fill_manual(values = c('purple','blue','green','orange'))
ggsave(here(paste0("figures/plots/kmk_seasonal_totkg.png")))

