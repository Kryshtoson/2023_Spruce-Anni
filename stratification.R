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
spe |> write_csv('species_data_long.csv')

spe |>
  distinct(rowname = PlotObservationID) |>
  left_join(read.table(r'(data\hlava.txt)', sep = '\t', header = T) |>
              as_tibble() |>
              select(rowname = `Relevé.number`, Lon = DEG_LON, Lat = DEG_LAT, cryptogams)) -> plot_coords
plot_cryptogams <- plot_coords[c('rowname', 'cryptogams')]
plot_coords <- plot_coords[1:3]

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
  geom_sf() +
  coord_sf(xlim = c(-10, 30),
           ylim = c(40, 70))

# tic()
# distance_plots <- plots_sf |> st_transform(3035) |> st_distance()
# colnames(distance_plots) <- plots_sf$rowname
# toc()
#
# step1 <- tibble(id = colnames(distance_plots)) |>
#   bind_cols(distance_plots |> as_tibble()) |>
#   pivot_longer(-1)
#
# step1 |>
#   mutate(value2 = as.numeric(value)) |>
#   filter(value2 < 5000) |>
#   group_by(id) |>
#   count() -> plots_around
# plots_around$n |> mean()
# hist(plots_around$n)
#
# plots_sf |>
#   left_join(plots_around |>
#               rename(rowname = id) |>
#               mutate(rowname = as.numeric(rowname))) |>
#   ggplot() +
#   geom_sf(data = world) +
#   geom_sf(aes(colour = n)) +
#   coord_sf(xlim = c(-10, 30),
#            ylim = c(40, 70))

#' ========================================================================
#' stratification
#' ========================================================================
plots_with_crypto <- plot_cryptogams$rowname[plot_cryptogams$cryptogams > 0]
coord_crypto <- coord[rownames(coord) %in% plots_with_crypto,]
spe_crypto <- spe |> filter(PlotObservationID %in% plots_with_crypto)

selection_crypto <- filtering(coord_crypto, spe_crypto,
                       dist.threshold = 10000,   #' 5 km
                       sim.threshold = 0.7,
                       sim.method = c("bray"),
                       remove = c("random"),
                       seed = 1234,
                       parallelize = TRUE)

selection_all <- filtering(coord, spe,
                       dist.threshold = 10000,   #' 5 km
                       sim.threshold = 0.7,
                       sim.method = c("bray"),
                       remove = c("random"),
                       seed = 1234,
                       parallelize = TRUE)


plots_sf |>
  dplyr::rename('plotID' = rowname) |>
  left_join(bind_rows(
  tibble(
    selection_run = 'with_crypto',
    plotID = as.numeric(rownames(selection_crypto))),
  tibble(
    selection_run = 'without_crypto',
    plotID = as.numeric(rownames(selection_all)))) |>
  mutate(value = T) |>
  pivot_wider(names_from = selection_run)) |>
  mutate(across(c('with_crypto', 'without_crypto'), function(x){ifelse(is.na(x), FALSE, x)})) |>
  write_sf('outputs//selection_output_10km_07Bray_4326.shp.gpkg')

read_sf('outputs//selection_output.gpkg') |>
  write_csv('outputs//selection_output_10km_07Bray.shp.csv')
