library(tidyverse)
library(oce)

# Load data
# NOTE: re-run collect_lunges.R when back from field to update DTAG deployment metadata
load("data/lunges_deployments.RData")
load("data/master_filtprey_data.RData")

# For each deployment:
# Species, location, prey, CEE start, day hours, night hours, day lunges, night lunges
# Deployment locations
deployment_locs <- readxl::read_xlsx("data/TAG GUIDE.xlsx", skip = 2) %>% 
  select(ID, lon = `Long_On    _`, lat = `Lat_On    _`) %>% 
  rbind(select(deployments_dtag, ID, lon = Long, lat = Lat)) %>% 
  mutate_at(vars(lon:lat), parse_double) %>% 
  drop_na() %>% 
  # fill in approximate lat/lon for two Antarctic deployments with missing location data
  rbind(tribble(~ID,           ~lon, ~lat,
                "mn170220-30", -65,  -65,
                "mn170322-30", -65,  -65))
deployments <- deployments_cats %>% 
  transmute(species, 
            ID, 
            tag = "CATS",
            tagon = tagon_prh, 
            tagoff = tagoff_prh, 
            sonar_exp = FALSE,
            cee_start = as.POSIXct(NA)) %>% 
  # Three deployments are missing tagoff (bw13_211a, bw15_076a, bw15_232a)
  rbind(deployments_dtag %>% 
          transmute(species, 
                    ID, 
                    tag = "DTAG",
                    tagon = `Tag on`, 
                    tagoff = `Tag off`,
                    sonar_exp,
                    cee_start)) %>% 
  # Add location
  left_join(select(deployment_locs, ID, lon, lat), by = "ID") %>% 
  # Add prey
  left_join(select(master_filtprey_data, ID, prey = Prey), by = "ID") %>% 
  # Remove a humpback bottom-feeding on inverts
  filter(ID != "mn151103-3") %>% 
  # Categorize prey as krill or fish
  mutate(prey = if_else(prey == "Krill", "Krill", "Fish")) %>% 
  # Bucket regions
  mutate(region = case_when(lat > 0 & lon < -90 ~ "NPac",
                            lat > 0 & lon > -90 ~ "NAtl",
                            lat < -60 & lon < 0 ~ "Ant",
                            lat < 0 & lon > -30 ~ "SAtl",
                            TRUE ~ NA_character_)) %>% 
  # Assign krill to two missing blue whales and Antarctic humpbacks
  mutate(prey = if_else(ID %in% c("bw15_076a", "bw15_232a", "mn10_133a", "mn10_139a", "mn10_144a", "mn10_151a", "mn10_155a"), 
                        "Krill",
                        prey),
         # Assign fish to two NPac humpbacks
         prey = if_else(ID %in% c("mn12_289a", "mn12_289b"),
                        "Fish",
                        prey))

# Pivot species and prey by region
group_by(deployments, region, tag, species, prey) %>% 
  summarize(N = n())  %>% 
  spread(region, N, fill = 0)

# Lunges for both tag types
lunge_tbl <- lunge_tbl_cats %>% 
  mutate(tag = "CATS") %>% 
  select(-tagon_guide,
         -tagoff_guide) %>% 
  rename(tagon = tagon_prh,
         tagoff = tagoff_prh) %>% 
  rbind(transmute(lunge_tbl_dtag,
                  species,
                  ID,
                  tagon = `Tag on`,
                  tagoff = `Tag off`,
                  lunge_dt,
                  tag = "DTAG"))

# Day and night
daynight_fun <- function(dt_UTC, lon, lat, night_thr = -18) {
  # This is a workaround because sunAngle can't handle NAs
  # Drop all NAs first, calculate sun angle, then join back onto original
  # data so the order of results is correct
  daynight <- tibble(dt_UTC, lon, lat) %>% 
    drop_na %>% 
    mutate(angle = sunAngle(dt_UTC, lon, lat)$altitude,
           daynight = factor(case_when(angle > 0 ~ 1,
                                       angle < night_thr ~ -1,
                                       TRUE ~ 0),
                             levels = c(1, -1, 0),
                             labels = c("Day", "Night", "Twilight")))
  result <- tibble(dt_UTC, lon, lat) %>% 
    left_join(daynight, by = c("dt_UTC", "lon", "lat")) 
  browser()
  result$daynight
}
daylight_fun <- function(start, end, lon, lat, n = 1e3) {
  pmap_dbl(list(start, end, lon, lat), 
           function(start, end, lon, lat) {
             if(is.na(lon) || is.na(lat)) {
               NA
             } else {
               angle <- sunAngle(seq(start, end, length.out = n), lon, lat)$altitude 
               as.numeric(end - start, units = "hours") * sum(angle > 0) / n
             }
           })
}
night_fun <- function(start, end, lon, lat, n = 1e3) {
  pmap_dbl(list(start, end, lon, lat), 
           function(start, end, lon, lat) {
             if(is.na(lon) || is.na(lat)) {
               NA
             } else {
               angle <- sunAngle(seq(start, end, length.out = n), lon, lat)$altitude 
               as.numeric(end - start, units = "hours") * sum(angle < -18) / n
             }
           })
}

# Lunges by day/night
lunge_daynight <- lunge_tbl %>%
  # Add locations to lunges
  left_join(deployment_locs, by = "ID") %>%
  # Figure out day/night and drop twilight
  mutate(daynight = daynight_fun(lunge_dt, lon, lat)) %>%
  filter(daynight != "Twilight") %>%
  group_by(ID, daynight) %>%
  summarize(N = n()) %>%
  ungroup %>%
  spread(daynight, N, fill = 0) %>%
  rename(day_lunges = Day,
         night_lunges = Night)
# Hours daylight
deploy_daynight <- ungroup(deployments_cats) %>% 
  mutate(tag = "CATS") %>% 
  select(-tagon_guide,
         -tagoff_guide) %>% 
  rename(tagon = tagon_prh,
         tagoff = tagoff_prh) %>% 
  rbind(transmute(ungroup(deployments_dtag),
                  species,
                  ID,
                  tagon = `Tag on`,
                  tagoff = `Tag off`,
                  tag = "DTAG")) %>%
  left_join(deployment_locs, by = "ID") %>%
  drop_na() %>% 
  mutate(daylight_hours = daylight_fun(tagon, tagoff, lon, lat),
         night_hours = night_fun(tagon, tagoff, lon, lat))

daily_rates <- left_join(deploy_daynight, lunge_daynight, by = "ID") %>% 
  mutate(day_rate = day_lunges / daylight_hours,
         night_rate = night_lunges / night_hours)

daily_rates %>% 
  filter(day_rate < 200, night_rate < 200) %>% 
  mutate(lunges_day = day_rate * 14.6 + night_rate * 5.5) %>% 
  ggplot(aes(species, lunges_day)) +
  geom_boxplot() +
  theme_minimal()

