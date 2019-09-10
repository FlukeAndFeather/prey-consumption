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

# Regex patterns for tag-type IDs
cats_ID <- "[a-z]{2}[0-9]{6}-[0-9]{1,2}[a-z]?"
dtag_ID <- "[a-z]{2}[0-9]{2}_[0-9]{1,3}[a-z]"
acou_ID <- "[a-z]{2}[0-9]{6}"

# Verified lunges
valid_deployments <- read_csv("data/lunge_rates_from_Paolo.csv") %>% 
  filter(lunge_quality %in% c("good", "good dives", "good_dives"), 
         sonar_exp == "none",
         species %in% c("ba", "bp", "bw", "mn"),
         total_lunges > 0) %>% 
  mutate(tag_type = case_when(str_detect(ID, cats_ID) ~ "CATS",
                              str_detect(ID, dtag_ID) ~ "DTAG",
                              str_detect(ID, acou_ID) ~ "ACOU"),
         species = recode(species, ba = "bb"))

# CATS tag guide
tag_guide <- readxl::read_xlsx("/Volumes/COPYCATSdat/CATS/TAG GUIDE.xlsx", skip = 2) %>% 
  mutate(tag_tz = sprintf("Etc/GMT%+d", `UTC         _`))

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
  tag_tz <- key$tag_tz
  data_paths <- find_cats_data(id)
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
  tag_tz <- key$tag_tz
  data_paths <- find_cats_data(id)
  tryCatch({
    prh <- readMat(data_paths$prh_path)
    tagonoff_prh <- dn_to_posix(range(prh$DN[prh$tagon == 1])) %>% 
      force_tz(tag_tz)
    tibble(tagon = tagonoff_prh[1],
           tagoff = tagonoff_prh[2])
  }, error = function(e) tibble(tagon = NA,
                                tagoff = NA))
}

# Deployment metadata. Starts with Paolo's lunge rate file. Filters down to
# good quality dives and the four rorqual species of interest. Removes sonar
# exposure deployments. Gets tag on/off times.
deployments_cats <- valid_deployments %>% 
  filter(tag_type == "CATS") %>% 
  left_join(select(tag_guide, ID, tag_tz), by = "ID") %>% 
  # One incorrect ID
  mutate(ID = recode(ID, `mn161106-36` = "mn161106-36b")) %>% 
  group_by(ID, tag_tz) %>% 
  group_map(tag_times_cats) %>% 
  ungroup

deployments_dtag <- read_xlsx("data/DTag data_fixesforMFC.xlsx", 
                              col_types = c("text", "date", "date", "numeric", 
                                            "numeric", "text", "numeric", 
                                            "numeric", "text", "text", 
                                            "text")) %>% 
  drop_na(ID, Lat, Long) %>%
  # Filter down to good quality lunges
  semi_join(valid_deployments %>% 
              filter(tag_type == "DTAG"),
            by = "ID") %>% 
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

