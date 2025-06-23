## Calculate Nicola Discharge
# Calculate August mean discharge following Warkentin et al. 2022
# June 23, 2025

library(tidyverse)
library(tidyhydat)


# id all stations
stations <- hy_stations()

# Look for stations with "Nicola" in the name
nicola_stations <- stations[
  grep("Nicola", stations$STATION_NAME, ignore.case = TRUE), 
] %>% 
  filter(HYD_STATUS == "ACTIVE") %>% 
  pull(STATION_NUMBER)


flow_data <- purrr::map(
  # no data for station 
  nicola_stations[c(1, 3, 4)], ~ hy_daily_flows(.x)
) %>% 
  bind_rows() %>% 
  mutate(
    year = lubridate::year(Date),
    day = lubridate::yday(Date),
    month = lubridate::month(Date) 
  ) %>% 
  filter(
    # year > 2010,
    month == "8"
  ) %>% 
  left_join(., 
            stations %>% 
              select(STATION_NUMBER, STATION_NAME),
            by = "STATION_NUMBER")

ggplot(flow_data) +
  geom_boxplot(
    aes(x = as.factor(day), y = Value)
  ) +
  facet_wrap(
    ~STATION_NAME
  )

flow_data %>% 
  group_by(
    year, STATION_NAME
  ) %>% 
  summarize(
    aug_flow = mean(Value)
  ) %>% 
  filter(
    STATION_NAME == "NICOLA RIVER AT OUTLET OF NICOLA LAKE"
  ) %>% 
  ggplot(.) +
  geom_point(aes(x = year, y = aug_flow)) 
