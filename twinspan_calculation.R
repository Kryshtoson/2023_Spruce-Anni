library(twinspan)
library(tidyverse)
#' ========================================================================
#' twinspan
#' ========================================================================

spe_ready <- read_csv('species_data_long.csv') |>
  semi_join(read_csv('stratification_selection.csv') |>
              filter(selected) |>
              select(PlotObservationID = rowname)) |>
  pivot_wider(names_from = Taxon_name, values_from = cover, values_fill = 0)

twin_out <- twinspan(spe_ready[, -1], cutlevels = c(3, 5, 25))

spe_ready |>
  distinct(PlotObservationID) |>
  left_join(read.table(r'(data\hlava.txt)', sep = '\t', header = T) |>
              as_tibble() |>
              select(PlotObservationID = `Relevé.number`, Lon = DEG_LON, Lat = DEG_LAT)) |>
  mutate(twin = cut(twin_out,1) |> factor()) -> plot_coords

plot_coords |>
  st_as_sf(coords = c('Lon', 'Lat'), crs = 4326) |>
  ggplot() +
  geom_sf(data = world) +
  geom_sf(aes(colour = twin)) +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))
