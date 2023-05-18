library(tictoc)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(twinspan)

world <- ne_countries(scale = "large", returnclass = "sf")

spe <- read.csv(r'(data\full_table.txt)', header = F) |>
  as_tibble() |>
  select(PlotObservationID = V1,
         Taxon_name = V2,
         cover = V4) |>
  filter(PlotObservationID != 32460,
         PlotObservationID != 77648) #' filtering Pyrenees

spe |>
  distinct(rowname = PlotObservationID) |>
  left_join(read.table(r'(data\hlava.txt)', sep = '\t', header = T) |>
              as_tibble() |>
              select(rowname = `Relevé.number`, Lon = DEG_LON, Lat = DEG_LAT)) -> plot_coords

coord <- plot_coords |>
  as.data.frame() |>
  column_to_rownames()

#' conversion to metric coordinate system
plots_sf <- plot_coords |> st_as_sf(coords = c('Lon', 'Lat'), crs = 4326)
plot_coords_laea <- plots_sf |>
  st_transform(3035) |>
  st_coordinates() |>
  as_tibble()
coord$Lat <- plot_coords_laea$X
coord$Lon <- plot_coords_laea$Y

#' extent of the original dataset
plots_sf |>
  ggplot() +
  geom_sf(data = world) +
  geom_sf(colour = 'red') +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))

tic()
distance_plots <- plots_sf |> st_transform(3035) |> st_distance()
colnames(distance_plots) <- plots_sf$rowname
toc()

step1 <- tibble(id = colnames(distance_plots)) |>
  bind_cols(distance_plots |> as_tibble()) |>
  pivot_longer(-1)

step1 |>
  mutate(value2 = as.numeric(value)) |>
  filter(value2 < 5000) |>
  group_by(id) |>
  count() -> plots_around
plots_around$n |> mean()
hist(plots_around$n)

plots_sf |>
  left_join(plots_around |>
              rename(rowname = id) |>
              mutate(rowname = as.numeric(rowname))) |>
  ggplot() +
  geom_sf(data = world) +
  geom_sf(aes(colour = n)) +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))

#' ========================================================================
#' stratification
#' ========================================================================

selection <- filtering(coord, spe,
                       dist.threshold = 10000,   #' 5 km
                       sim.threshold = 0.7,
                       sim.method = c("bray"),
                       remove = c("random"),
                       seed = 1234,
                       parallelize = TRUE)

plots_sf |>
  left_join(tibble(rowname = as.numeric(rownames(selection)),
                   selected = T)) |>
  mutate(selected = replace_na(selected, F)) -> selection_stats

selection_stats$selected |> table()
selection_stats |> write_sf('spatial_outputs\\stratified_10km_07Bray_4326.shp')

selection_stats |> ggplot() +
  geom_sf(data = world) +
  geom_sf(aes(colour = selected)) +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))

selection_stats |> write_csv('stratification_selection.csv')

#' ========================================================================
#' twinspan
#' ========================================================================

spe_ready <- spe |>
  semi_join(selection_stats |>
    as_tibble() |>
              filter(selected)|>
              select(PlotObservationID = rowname)) |>
  pivot_wider(names_from = Taxon_name, values_from = cover, values_fill = 0)

twin_out <- twinspan(spe_ready[, -1], cutlevels = c(3, 5, 25))

selection_stats |>
  filter(selected) |>
  mutate(twin = cut(twin_out, 2) |> factor()) |>
  ggplot() +
  geom_sf(data = world) +
  geom_sf(aes(colour = twin)) +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))
