library(lubridate)
library(lutz)
library(oce)
library(R.matlab)
library(readxl)
library(tidyverse)
library(zeallot)

# Convert MATLAB DN's to POSIXct
dn_to_posix <- function(dn) {
  as.POSIXct((dn - 719529) * 86400, origin = "1970-01-01", tz = "UTC")
}

# What are the day/night feeding rates?
# Look up PRH & lunge file paths on the CATS drive
# This will look different on Windows
find_cats_data <- function(id) {
  # Look for CATS
  cats_dir <- dir("/Volumes/COPYCATSdat/CATS/tag_data", 
                  pattern = id, 
                  full.names = TRUE)
  # There's one deployment where the lunges mat file isn't in a separate folder
  lunge_dir = if (id == "bw140806-2") {
    cats_dir
  } else {
    file.path(cats_dir, "lunges")
  }
  list(prh_path = dir(cats_dir,
                      pattern = "10Hzprh",
                      full.names = TRUE),
       lunge_path = dir(lunge_dir,
                        pattern = "lunges.*mat",
                        full.names = TRUE))
}
find_dtag_data <- function(id) {
  # Look for dtag
  dtag_dir <- "/Volumes/COPYCATSdat/DTAG"
  list(prh_path = dir(file.path(dtag_dir, "prh"),
                      pattern = sprintf("%s.*", id),
                      full.names = TRUE),
       lunge_path = dir(file.path(dtag_dir, "lunges"),
                        pattern = sprintf("%s.*", id),
                        full.names = TRUE))
}

# Find times of lunges. To be used with group_map. Takes an ID, finds the
# associated PRH and lunge files, converts times to POSIX, and returns a tibble
# with lunge times.d
lunge_times_cats <- function(data, key) {
  id <- key$ID
  data_paths <- find_cats_data(id)
  tag_guide <- readxl::read_xlsx("/Volumes/COPYCATSdat/CATS/TAG GUIDE.xlsx", skip = 2)
  tag_tz <- sprintf("Etc/GMT%+d", -filter(tag_guide, ID == id)$`UTC         _`)
  tryCatch({
    c(prh, lunges) %<-% map(data_paths, readMat)
    tibble(lunge_dt = dn_to_posix(prh$DN[lunges$LungeI]) %>% 
             force_tz(tag_tz))
  },
  error = function(e) tibble(lunge_dt = NA))
}
lunge_times_dtag <- function(data, key) {
  id <- key$ID
  tagon <- key$`Tag on`
  data_paths <- find_dtag_data(id)
  tryCatch({
    c(prh, lunges) %<-% map(data_paths, readMat)
    fs <- prh$fs[1]
    tibble(lunge_dt = (lunges$LungeI - 1) * fs + tagon)
  },
  error = function(e) tibble(lunge_dt = NA))
}

# Finds tag on and tag off times. To be used with group_map. Takes an ID, loads
# the PRH, finds the tag on/off times, and converts to POSIX. Returns a tibble
# with two columns (tagon and tagoff).
tag_times_cats <- function(data, key) {
  id <- key$ID
  data_paths <- find_cats_data(id)
  tryCatch({
    # Read time zone from tag guide
    tag_guide <- readxl::read_xlsx("/Volumes/COPYCATSdat/CATS/TAG GUIDE.xlsx", skip = 2)
    tag_tz <- sprintf("Etc/GMT%+d", -filter(tag_guide, ID == id)$`UTC         _`)
    prh <- readMat(data_paths$prh_path)
    tagonoff_prh <- dn_to_posix(range(prh$DN[prh$tagon == 1])) %>% 
      force_tz(tag_tz)
    # Tag_On seems to read correctly, but Tag_Off comes across as a serial
    tagonoff_guide <- filter(tag_guide, ID == id)$Tag_On %>% 
      force_tz(tag_tz)
    tagonoff_guide[2] <- filter(tag_guide, ID == id)$Tag_Off %>% 
      openxlsx::convertToDateTime(tz = tag_tz)
    tibble(tagon_prh = tagonoff_prh[1],
           tagoff_prh = tagonoff_prh[2],
           tagon_guide = tagonoff_guide[1],
           tagoff_guide = tagonoff_guide[2])
  }, error = function(e) tibble(tagon_prh = NA,
                                tagoff_prh = NA,
                                tagon_guide = NA,
                                tagoff_guide = NA))
}

# Deployment metadata. Starts with Paolo's lunge rate file. Filters down to
# good quality dives and the four rorqual species of interest. Removes sonar
# exposure deployments. Gets tag on/off times.
deployments_cats <- read_csv("data/lunge_rates_from_Paolo.csv") %>% 
  filter(lunge_quality %in% c("good", "good dives", "good_dives"), 
         sonar_exp == "none",
         species %in% c("ba", "bp", "bw", "mn")) %>% 
  # One incorrect ID
  mutate(ID = if_else(ID == "mn161106-36",
                      "mn161106-36b",
                      ID)) %>% 
  group_by(species, ID) %>% 
  group_map(tag_times_cats) %>% 
  ungroup %>% 
  filter(str_length(ID) >= 11)

deployments_dtag <- read_xlsx("data/DTag data_fixesforMFC.xlsx", 
                              col_types = c("text", "date", "date", "numeric", 
                                            "numeric", "text", "numeric", 
                                            "numeric", "text", "text", 
                                            "text")) %>% 
  drop_na(ID, Lat, Long) %>%
  mutate(species = substr(ID, 1, 2),
         deploy_tz = tz_lookup_coords(lat = Lat, lon = Long, method = "accurate"),
         tagon = force_tzs(`Tag on`, deploy_tz),
         tagoff = force_tzs(`Tag off`, deploy_tz),
         cee_start = ISOdatetime(year(`Tag on`), 
                                 month(`Tag on`), 
                                 day(`Tag on`), 
                                 ifelse(nchar(`Time CEE start`) == 4, 
                                        substr(`Time CEE start`, 1, 2), 
                                        substr(`Time CEE start`, 1, 1)), 
                                 ifelse(nchar(`Time CEE start`) == 4, 
                                        substr(`Time CEE start`, 3, 4), 
                                        substr(`Time CEE start`, 2, 3)), 
                                 0) %>% 
           force_tzs(deploy_tz))

# Looks up the lunge times for all the deployments.
lunge_tbl_cats <- deployments_cats %>% 
  group_by_all %>% 
  group_map(lunge_times_cats) %>% 
  ungroup
lunge_tbl_dtag <- deployments_dtag %>% 
  group_by_all %>% 
  group_map(lunge_times_dtag) %>% 
  ungroup

save(deployments_cats, deployments_dtag, lunge_tbl_cats, lunge_tbl_dtag, 
     file = "data/lunges_deployments.RData")

