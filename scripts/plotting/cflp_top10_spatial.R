
rm(list=ls())
gc()
library(cmocean)
library(dplyr)
library(lubridate)
library(sf)
library(terra)
library(viridisLite)
library(ggplot2)

setwd("C:/Users/brendan.turley/Documents/data/shapefiles/GSHHS_shp/i")
world <- vect('GSHHS_i_L1.shp') |> 
  crop(ext(-100,-80,23,32))

setwd("C:/Users/brendan.turley/Documents/data/shapefiles/king_mackerel")
kmk <- vect('king_mackerel_po.shp')
# setwd("C:/Users/brendan.turley/Documents/data/shapefiles/GOM_2500ft")
# gom <- vect('GOM_2500ft.shp')

setwd("C:/Users/brendan.turley/Documents/data/shapefiles/cflp_statgrid")
sz_shp <- vect('CFLP_StatGrid_2013_v20140210.shp') |>
  st_as_sf()
sz_shp$AREA_FISHED <- sz_shp$SZ_ID

#### read data and subset ####--------------------------------------------------
gom_st <- c('FL', 'AL', 'MS', 'LA', 'TX')

setwd("C:/Users/brendan.turley/Documents/CMP/data/cflp")
cflp <- readRDS('CFLPblake.rds')
yr <- 2013 # 2019
cflp <- subset(cflp, LAND_YEAR>yr & CATCH_TYPE == 'CATCH') |>
  subset(COMMON_NAME=='MACKERELS, KING AND CERO') |>
  subset(REGION == 'GOM' & is.element(ST_ABRV, gom_st))  |>
  subset(FLAG_MULTIGEAR==0 & FLAG_MULTIAREA==0 & !is.na(AREA_FISHED))
gc()


# find top 10 areas by total landings per year
tot_lbs_yr_shp <- aggregate(TOTAL_WHOLE_POUNDS ~ LAND_YEAR + AREA_FISHED,
                            data = cflp, sum, na.rm = T)

kmk_areas <-  merge(tot_lbs_yr_shp,
                    sz_shp,
                    by = c('AREA_FISHED')) |>
  st_as_sf()

area_top10 <- kmk_areas |>
  filter(LAND_YEAR >= 2020) |>
  group_by(LAND_YEAR) |> 
slice_max(order_by = TOTAL_WHOLE_POUNDS, n = 10)


# for(i in 2020:2024){
#   vect(kmk_areas) |> plot()
#   subset(kmk_areas, LAND_YEAR==i) |>
#     # arrange(desc(TOTAL_WHOLE_POUNDS)) |>
#     slice_max(order_by = TOTAL_WHOLE_POUNDS, n = 10) |> 
#     vect() |> plot(y='TOTAL_WHOLE_POUNDS', add = T)
#   mtext(i)
# }

### find top 10 areas by number of trip per year
trips_yr_shp <- aggregate(SCHEDULE_NUMBER ~ LAND_YEAR + AREA_FISHED,
                            data = cflp, function(x) length(unique(x)))

trips_areas <-  merge(trips_yr_shp,
                    sz_shp,
                    by = c('AREA_FISHED')) |>
  st_as_sf()

trips_top10 <- trips_areas |>
  filter(LAND_YEAR >= 2020) |>
  group_by(LAND_YEAR) |> 
  slice_max(order_by = SCHEDULE_NUMBER, n = 10)

# for(i in 2020:2024){
#   vect(trips_areas) |> plot()
#   subset(trips_areas, LAND_YEAR==i) |>
#     # arrange(desc(TOTAL_WHOLE_POUNDS)) |>
#     slice_max(order_by = SCHEDULE_NUMBER, n = 10) |>
#     vect() |> plot(y='SCHEDULE_NUMBER', add = T)
#   mtext(i)
# }

### aggregate exploratory plots
tot_lbs_shp <- aggregate(TOTAL_WHOLE_POUNDS ~ AREA_FISHED,
                            data = cflp, sum, na.rm = T) |> 
  merge(sz_shp, by = c('AREA_FISHED')) |>
  st_as_sf()

mean_per_year <- cflp |>
  group_by(AREA_FISHED, LAND_YEAR) |>
  summarize(yearly_mean = mean(TOTAL_WHOLE_POUNDS, na.rm = TRUE), .groups = "drop") |>
  group_by(AREA_FISHED) |>
  summarize(overall_mean = mean(yearly_mean, na.rm = TRUE), .groups = "drop") |> 
  merge(sz_shp, by = c('AREA_FISHED')) |>
  st_as_sf()

trips_shp <- aggregate(SCHEDULE_NUMBER ~ AREA_FISHED,
                          data = cflp, function(x) length(unique(x))) |>
  merge(sz_shp, by = c('AREA_FISHED')) |>
  st_as_sf()

trips_per_year <- cflp |>
  group_by(AREA_FISHED, LAND_YEAR) |>
  summarize(yearly_trips = length(unique(SCHEDULE_NUMBER)), .groups = "drop") |>
  group_by(AREA_FISHED) |>
  summarize(overall_trips = mean(yearly_trips, na.rm = TRUE), .groups = "drop") |> 
  merge(sz_shp, by = c('AREA_FISHED')) |>
  st_as_sf()

ggplot(data = tot_lbs_shp) +
  geom_sf(aes(fill = TOTAL_WHOLE_POUNDS), color = "white", size = 0.2) +
  # facet_wrap(~ LAND_YEAR, ncol = 2) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24,30.5), xlim = c(-97, -81)) +
  labs(fill = 'Total landings (lbs)') +
  theme_bw()

ggplot(data = mean_per_year) +
  geom_sf(aes(fill = overall_mean), color = "white", size = 0.2) +
  # facet_wrap(~ LAND_YEAR, ncol = 2) +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24,30.5), xlim = c(-97, -81)) +
  labs(fill = 'Mean landings (lbs)') +
  theme_bw()

ggplot(data = trips_shp) +
  geom_sf(aes(fill = SCHEDULE_NUMBER), color = "white", size = 0.2) +
  # facet_wrap(~ LAND_YEAR, ncol = 2) +
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24,30.5), xlim = c(-97, -81)) +
  labs(fill = 'Number of Trips') +
  theme_bw()

ggplot(data = trips_per_year) +
  geom_sf(aes(fill = overall_trips), color = "white", size = 0.2) +
  # facet_wrap(~ LAND_YEAR, ncol = 2) +
  scale_fill_viridis_c(option = "rocket", direction = -1) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24,30.5), xlim = c(-97, -81)) +
  labs(fill = 'Mean Number of Trips') +
  theme_bw()


### plots for ESP
ggplot(data = area_top10) +
  geom_sf(aes(fill = TOTAL_WHOLE_POUNDS), color = "white", size = 0.2) +
  facet_wrap(~ LAND_YEAR, ncol = 2) + 
  scale_fill_viridis_c(option = "mako", direction = -1, labels = label_comma()) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24.5,30.5), xlim = c(-97, -81.5)) +
  labs(fill = 'Landings (lbs)') +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.65, .15),
        legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 10))
ggsave('top10_landings.png', width = 7, height = 5, units = 'in',
       path = "~/R_projects/King-Mackerel-ESP/figures/plots")

ggplot(data = trips_top10) +
  geom_sf(aes(fill = SCHEDULE_NUMBER), color = "white", size = 0.2) +
  facet_wrap(~ LAND_YEAR, ncol = 2) + 
  scale_fill_viridis_c(option = "rocket", direction = -1, labels = label_comma()) +
  geom_sf(data = world |> st_as_sf()) +
  coord_sf(ylim = c(24.5,30.5), xlim = c(-97, -81.5)) +
  labs(fill = 'Trips') +
  theme_bw() +
  theme(legend.position = "inside", legend.position.inside = c(0.65, .15),
        legend.key.size = unit(0.5, "cm"), legend.title = element_text(size = 10))
ggsave('top10_trips.png', width = 7, height = 5, units = 'in',
       path = "~/R_projects/King-Mackerel-ESP/figures/plots")


### pull out handlines only
# table(cflp$GEAR)
gear_keep <- c('H', 'E', 'TR')
# gear_keep <- c('TR')
cflp_hl <- subset(cflp , is.element(cflp$GEAR, gear_keep)) |>
  subset(FLAG_MULTIGEAR==0 & FLAG_MULTIAREA==0)