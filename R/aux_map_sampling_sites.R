# this script plots a map of the study sites

# INPUTS
# - CBMS transects and visits (all_visits.txt, transectos_CBMS_filtrados.txt)
# - UBMS transects and visits (all_visits.txt, ubms_sites.csv)
# - land cover raster (ESA_WorldCover_10m_2021_V200_N39E000_Map.tif)
# - built surface raster (GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R5_C19.tif)

# OUTPUTS
# - figure (bcn_urb_map)

# NOTE: This script will not work directly, as it requires data that is 
# not available in the open repository (the folder "map_data"). See the section on "code and data availability" 
# of the main text for details.

# NOTE: 2 this code is adapted from a script by Pau Colom. It can be tweaked
# to plot more complex figures, e.g. with land uses or built-up surface in the background
# but these are not needed here, as the objective is to show the spatial arrangement of
# sampling points

# -------------------------------------------------------------------------
# Clean env and load packages ---
library(tidyverse)
library(sf)
library(raster)
library(osmdata)
library(viridis)
library(ggspatial)
library(terra)
library(patchwork)
# -------------------------------------------------------------------------
# Load data ---

# NOTE: this script will not work out of the box as it requires data that 
# is not uploaded to ZENODO . See the README for details on files uploaded and permissions

# The path to the data used is set here:
data_path <- "~/trabajo/projects/BCN/butterfly_trait_distributions/"

# rasters
lc_raster <- rast(paste0(data_path,"data/map_data/ESA_WorldCover_10m_2021_V200_N39E000_Map.tif"))
ghsl <- raster(paste0(data_path,"data/map_data/GHS_BUILT_S_E2020_GLOBE_R2023A_4326_3ss_V1_0_R5_C19.tif"))

# sites
all_sites <- read.csv2(paste0(data_path,"data/map_data/all_visits.txt")) %>% dplyr::select(SITE_ID,SCHEME) %>% 
  unique

cbms_coords <- read.csv2(paste0(data_path,"data/map_data/transectos_CBMS_filtrados.csv"),sep=",") %>%
  dplyr::select(SITE_ID,lat,lng)
names(cbms_coords) <- c("SITE_ID","latitude","longitude")
cbms_coords$SITE_ID <- as.character(cbms_coords$SITE_ID)

ubms_coords <- read.csv2(paste0(data_path,"data/map_data/ubms_sites.csv")) %>%
  dplyr::select(transect_id, transect_latitude, transect_longitude)
names(ubms_coords) <- c("SITE_ID","latitude","longitude")
ubms_coords <- ubms_coords[grepl("_A", ubms_coords$SITE_ID),]
ubms_coords$SITE_ID <- substr(ubms_coords$SITE_ID, 1,nchar(ubms_coords$SITE_ID)-2)

cbms_sites <- as_tibble(all_sites) %>%
  left_join(cbms_coords) %>% drop_na()
ubms_sites <- as_tibble(all_sites) %>%
  left_join(ubms_coords) %>% drop_na()

all_sites <- bind_rows(cbms_sites,ubms_sites)

# --- Convert to sf and merge ---
# some of this is likely repetitive as it comes from reused code

ubms_sf <- st_as_sf(ubms_sites,
                    coords = c("longitude", "latitude"),
                    crs = 4326, remove = FALSE)

ubms_df <- ubms_sf %>%
  st_drop_geometry() %>%
  mutate(source = "uBMS") %>%
  dplyr::select(SITE_ID, longitude, latitude, source)

cbms_df <- cbms_sites %>%
  mutate(source = "CBMS") %>%
  dplyr::select(SITE_ID, longitude, latitude, source)

merged_transects <- bind_rows(cbms_df, ubms_df)

m_coord_clean <- merged_transects %>%
  filter(!is.na(longitude), !is.na(latitude))

m_coord_sf_wgs <- st_as_sf(m_coord_clean,
                       coords = c("longitude", "latitude"),
                       crs = 4326, remove = FALSE)

# --- Get Barcelona boundary ---

bcn_boundary <- opq("Barcelona") %>%
  add_osm_feature(key = "name", value = "Barcelona") %>%
  osmdata_sf() %>%
  {.$osm_multipolygons} %>%
  .[1, ]

# --- Extract built-up values + inside/outside classification ---
bcn_sites <- st_transform(m_coord_sf_wgs, crs = crs(ghsl))
bcn_sites$urban_value <- raster::extract(ghsl, bcn_sites)

bcn_sites$context <- ifelse(
  as.logical(st_within(bcn_sites, bcn_boundary, sparse = FALSE)[,1]),
  "inside", "outside"
)
bcn_sites$context <- factor(bcn_sites$context, levels = c("inside","outside"))

# --- Extract urban values for sites ---
bcn_sites <- st_transform(m_coord_sf_wgs, crs = crs(ghsl))
bcn_sites$urban_value <- raster::extract(ghsl, bcn_sites)

# --- Buffer around Barcelona boundary ---
bcn_boundary_utm <- st_transform(bcn_boundary, 25831) # ETRS89 / UTM zone 31N per Catalunya
bcn_buffer <- st_buffer(bcn_boundary_utm, dist = 60000) # 60 km buffer
bcn_buffer_wgs <- st_transform(bcn_buffer, crs = crs(ghsl))
bcn_buffer_sp <- as_Spatial(bcn_buffer_wgs)

# --- Crop/mask raster ---
ghsl_crop <- crop(ghsl, bcn_buffer_sp)
ghsl_mask <- mask(ghsl_crop, bcn_buffer_sp)

ghsl_df_bcn <- as.data.frame(rasterToPoints(ghsl_mask))
names(ghsl_df_bcn)[3] <- "layer"
ghsl_df_bcn$log_layer <- log1p(ghsl_df_bcn$layer)

# --- Inside vs outside classification ---
inside_mat <- st_within(bcn_sites, bcn_boundary, sparse = FALSE)
inside_flag <- as.logical(inside_mat[,1])

# bcn_sites$context <- ifelse(inside_flag, "inside", "outside")
# bcn_sites$context <- factor(bcn_sites$context, levels = c("inside", "outside"))

bcn_sites$context <- ifelse(inside_flag, "filtered", "regional")
bcn_sites$context <- factor(bcn_sites$context, levels = c("filtered", "regional"))

# --- Filter sites within buffer ---
bcn_sites$in_buffer <- as.logical(st_within(bcn_sites, st_as_sf(bcn_buffer_wgs), sparse = FALSE)[,1])
sites_buffer_bcn <- bcn_sites[bcn_sites$in_buffer, ]

# --- Reprojectar transectes al CRS del raster ---
bcn_sites <- st_transform(m_coord_sf_wgs, crs = crs(lc_raster))

# Extreure valors land cover per cada punt
bcn_sites$lc_value <- terra::extract(lc_raster, vect(bcn_sites))[,2]

# --- Crop i mask al buffer ---
lc_crop <- terra::crop(lc_raster, vect(bcn_buffer_wgs))
lc_mask <- terra::mask(lc_crop, vect(bcn_buffer_wgs))

lc_mask_filled <- lc_mask
lc_mask_filled[is.na(lc_mask_filled)] <- 80

# --- Reduir resolució per fer-ho més lleuger (250 m) ---
lc_raster_low <- aggregate(lc_mask_filled, fact = 25, fun = modal)
lc_df_bcn <- as.data.frame(lc_raster_low, xy = TRUE)
names(lc_df_bcn)[3] <- "class"

lc_df_bcn$class <- factor(
  lc_df_bcn$class,
  levels = c("50", "10", "20", "30", "60", "40", "90", "80")  # ordre desitjat
)

water_df_bcn <- lc_df_bcn %>% # Run BCN land cover map
  dplyr::filter(class == 80)

# -------------------------------------------------------------------------
# Plot

map_urb_bcn <- ggplot() +
  # geom_raster(data = ghsl_df_bcn, aes(x = x, y = y, fill = log_layer), alpha = 0.3) +
  geom_raster(data = ghsl_df_bcn, aes(x = x, y = y), fill = "darkolivegreen3", alpha = 0.3) +
  geom_sf(data = bcn_boundary, fill = NA, color = "grey30", size = 1) +
  geom_raster(data = water_df_bcn, aes(x = x, y = y), fill = "white", alpha = 1) +  # <-- water layer
  geom_sf(data = bcn_buffer_wgs, fill = NA, size = .7) +
  geom_sf(data = sites_buffer_bcn, aes(color = context), size = 1.5, alpha = 0.9) +
  scale_fill_viridis(option = "D", na.value = NA, name = "Built-up surface") +
  scale_color_manual(name = "Community", values = c("filtered" = "darkred", "regional" = "black")) +
  scale_x_continuous(
    breaks = seq(1.5, 3, by = 0.5)  
  ) +
  scale_y_continuous(
    breaks = seq(41, 42, by = 0.5),limits = c(41.18,42)
  )+
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_family = "Garamond", base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(family = "Garamond", size = 14),
    axis.title = element_text(family = "Garamond", size = 18),
    legend.position=c(.13, .88),
    legend.text = element_text(family = "Garamond", size = 12),
    legend.title = element_text(family = "Garamond", size = 14)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  annotation_scale(
    location = "tr", 
    width_hint = 0.25,
    text_family = "Garamond",  
    text_cex = 1            
  )

# map_urb_bcn

# -------------------------------------------------------------------------
# zoom over the city
map_city <- ggplot() +
  geom_sf(data = bcn_boundary, fill = NA, color = "grey30", linewidth = .4) +
  # geom_raster(data = water_df_bcn, aes(x = x, y = y), fill = "white", alpha = 1) +  # <-- water layer
  # geom_sf(data = bcn_buffer_wgs, fill = NA, size = .7) +
  geom_sf(data = subset(sites_buffer_bcn, context == "filtered"), color = "darkred", size = 1.8, alpha = 0.9) +
  lims(x = c(2.05,2.24), y = c(41.32,41.47)) +
  scale_x_continuous(
    # breaks = seq(1.5, 3, by = 0.5)
    breaks = NULL
  ) +
  scale_y_continuous(
    # breaks = seq(41, 42, by = 0.5),limits = c(41.18,42)
    breaks = NULL
  )+
  # labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_family = "Garamond", base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(family = "Garamond", size = 14),
    axis.title = element_text(family = "Garamond", size = 18),
    legend.position=c(.15, .9),
    legend.text = element_text(family = "Garamond", size = 14),
    legend.title = element_text(family = "Garamond", size = 16)
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  # annotation_scale(
  #   location = "tr", 
  #   width_hint = 0.25,
  #   text_family = "Garamond",  
  #   text_cex = 1            
  # ) + 
  NULL

# map_city

map_full <- map_urb_bcn + map_city

# -------------------------------------------------------------------------
# store

# ggsave("results/images/bcn_urb_map.pdf",
#        map_full,
#        device = cairo_pdf,
#        width = 10,
#        height = 5,
#        dpi = 300,
#        bg = "white")

