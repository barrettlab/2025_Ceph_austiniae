# Mapping Cephalanthera distibutions

# Load required libraries

# May need to install Polychrome this way
# install.packages("https://cran.r-project.org/src/contrib/Archive/Polychrome/Polychrome_1.4.0.tar.gz", repos = NULL, type = "source")

library(rgbif)
library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(Polychrome)
library(ggsci)

# --------------------------
# 1. Search GBIF for genus Cephalanthera
name_data <- name_backbone(name = "Cephalanthera", rank = "genus")
taxon_key <- name_data$usageKey

# 2. Download occurrence data (up to 20,000 points)
ceph_data <- occ_search(taxonKey = taxon_key,
                        hasCoordinate = TRUE,
                        limit = 20000)
						
# 3. Convert to dataframe
df <- ceph_data$data

# 4. Clean and filter the data
df_clean <- df %>%
  select(species, decimalLongitude, decimalLatitude, countryCode) %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude), !is.na(species))

# --------------------------
# 5. Convert to spatial object
points_sf <- st_as_sf(df_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# 6. Load world basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# --------------------------
# 7. Choose color palette
n_species <- length(unique(df_clean$species))

# if (n_species <= 9) {
#   # If 9 or fewer species, use bright ggplot palette
#   palette_choice <- scale_color_npg(name = "Species")
# } else if (n_species <= 20) {
#   # If 10–20 species, use a colorful scientific palette
#   palette_choice <- scale_color_d3(name = "Species")
# } else {
#   # If >20 species, use highly distinct Polychrome
#   palette_choice <- scale_color_manual(values = palette36.colors(n_species), name = "Species")
# }

# 8. Set up color palette correctly
species_names <- unique(df_clean_thinned$species)

palette_colors <- setNames(palette36.colors(length(species_names)), species_names)

palette_choice <- scale_color_manual(
  values = palette_colors,
  name = "Species"
)

# --------------------------
# 8. Plot!
p <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", size = 0.2) +
  geom_sf(data = points_sf, aes(color = species), size = 2, alpha = 0.8) +
  palette_choice +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray90"),
    axis.title = element_blank()
  ) +
  labs(
    title = "Global Distribution of *Cephalanthera* Species",
    subtitle = "Data from GBIF",
    caption = "Map created in R with rgbif + ggplot2"
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

# Print map
print(p)

# --------------------------
# 9. Save high-resolution figure
ggsave("Cephalanthera_species_map.png", p, width = 12, height = 8, dpi = 600)