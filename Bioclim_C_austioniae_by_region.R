# Data from GBIF via rgbif (Chamberlain et al. 2023. https://doi.org/10.21105/joss.00753)
# WorldClim bioclimatic data via geodata (Hijmans et al. 2023. https://rspatial.org/)
# Spatial raster handling via terra (Hijmans 2023. https://rspatial.org/)
# Spatial vector handling via sf (Pebesma 2018. https://doi.org/10.32614/RJ-2018-009)
# Data manipulation via dplyr (Wickham et al. 2023. https://dplyr.tidyverse.org/)
# Plotting via ggplot2 (Wickham 2016. https://ggplot2.tidyverse.org/)
# Biplot annotation via ggrepel (Slowikowski 2021. https://github.com/slowkow/ggrepel)
# PCA and matrix operations via base R stats (R Core Team 2024. https://www.r-project.org/)
# Analysis and script supported by OpenAI's ChatGPT-4 (2024. https://openai.com/chatgpt)



# ---- 1. Load packages ----
install.packages(c("rgbif", "geodata", "terra", "sf", "dplyr", "ggplot2"))
library(rgbif)
library(geodata)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ---- 2. Download GBIF data ----
occ_data <- occ_search(
  scientificName = "Cephalanthera austiniae",
  limit = 10000,
  hasCoordinate = TRUE
)

occ_df <- occ_data$data %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  dplyr::select(decimalLatitude, decimalLongitude, stateProvince, countryCode)

occ_df <- dplyr::filter(occ_data$data, !is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  dplyr::select(decimalLatitude, decimalLongitude, stateProvince, countryCode)

# ---- 3. Convert to spatial points ----
occ_sf <- st_as_sf(occ_df, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

# ---- 4. Download WorldClim bioclim data ----
wc_all <- geodata::worldclim_global(var = "bio", res = 10, path = tempdir())

# ---- 4b. Select 5 uncorrelated bioclim variables ----
# Extract all 19 BIO variables temporarily
wc19 <- wc_all

# Extract raster values at occurrence points
bio19_vals <- terra::extract(wc19, vect(occ_sf))

# Remove non-BIO columns & rows with missing data
bio_only <- bio19_vals[, grep("^bio", names(bio19_vals), ignore.case = TRUE)]
bio_only <- bio_only[complete.cases(bio_only), ]

# Compute correlation matrix
cor_mat <- cor(bio_only, use = "pairwise.complete.obs")

# Identify and remove highly correlated variables
# (change cutoff from 0.7 → 0.8 or 0.9 depending on your needs)
library(caret)
high_corr <- findCorrelation(cor_mat, cutoff = 0.70, exact = TRUE)

# Variables to keep
vars_keep <- names(bio_only)[-high_corr]

# If more than 5 remain, keep the first 5 for PCA
vars_final <- vars_keep[1:5]

cat("Selected variables with low correlation:\n")
print(vars_final)

# Subset the WC stack to selected variables
wc <- wc19[[vars_final]]
names(wc) <- toupper(vars_final)   # Make names BIO1, BIO2, etc.


names(wc) <- paste0("BIO", bioclim_vars)

# ---- 5. Extract bioclim values at occurrence points ----
occ_clim <- terra::extract(wc, vect(occ_sf)) %>%
  bind_cols(occ_df)

# ---- 6. Convert to data frame and clean ----
occ_df2 <- as.data.frame(occ_clim)

# Remove rows with missing bioclim data
bio_cols <- grep("^BIO", names(occ_df2), value = TRUE)
pca_input <- scale(occ_df2[, bio_cols])
valid_rows <- which(complete.cases(pca_input) & rowSums(is.finite(pca_input)) == ncol(pca_input))
pca_input_clean <- pca_input[valid_rows, ]
occ_df_clean <- occ_df2[valid_rows, ]

# ---- 7. Assign subregions ----
assign_subregion <- function(lat, lon, state) {
  if (is.na(state)) {
    return(NA)
  } else if (state == "California") {
    if (lon < -121.5) return("CA-Coast/Cascades") else return("CA-Sierra")
  } else if (state == "Oregon") {
    if (lon < -120.5) return("OR-West") else return("OR-East")
  } else if (state == "Washington") {
    if (lon < -121.5) return("WA-West") else return("WA-East")
  } else if (state == "Idaho") {
    return("ID")
  } else if (state == "British Columbia") {
    return("BC")
  } else {
    return("Other")
  }
}


occ_df_clean$Subregion <- mapply(assign_subregion, 
                                  occ_df_clean$decimalLatitude, 
                                  occ_df_clean$decimalLongitude, 
                                  occ_df_clean$stateProvince)

# ---- 8. Run PCA ----
pca_res <- prcomp(pca_input_clean)

# ---- 9. Plot PCA ----
scores <- as.data.frame(pca_res$x)
scores$Subregion <- occ_df_clean$Subregion

ggplot(scores, aes(x = PC1, y = PC2, color = Subregion)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA of Bioclim Variables for Cephalanthera austiniae",
       x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)")) +
  scale_color_brewer(palette = "Dark2")


# ---- 10. Plot to include biplot and convex hulls

library(ggrepel)

# ---- 1. Prepare biplot arrows (scaled loadings) ----
loadings <- pca_res$rotation[, 1:2]
loadings_df <- as.data.frame(loadings)
loadings_df$Variable <- rownames(loadings_df)

# Scale arrows for display (adjust multiplier if needed)
arrow_scale <- 3  # adjust this factor as needed
loadings_df <- loadings_df %>%
    mutate(PC1 = PC1 * arrow_scale,
           PC2 = PC2 * arrow_scale)

# ---- 2. Prepare convex hulls for subregions ----
hull_df <- scores %>%
    group_by(Subregion) %>%
    slice(chull(PC1, PC2))

# ---- 3. Plot PCA with biplot and convex hulls ----
ggplot(scores, aes(x = PC1, y = PC2, color = Subregion)) +
    # Convex hulls
    geom_polygon(data = hull_df, aes(fill = Subregion), color = NA, alpha = 0.1) +
    # Points
    geom_point(size = 2, alpha = 0.95) +
    # Biplot arrows
    geom_segment(data = loadings_df,
                 aes(x = 0, y = 0, xend = PC1, yend = PC2),
                 arrow = arrow(length = unit(0.25, "cm")),
                 inherit.aes = FALSE, color = "black", linewidth = 0.8) +
    # Variable labels
    geom_text_repel(data = loadings_df,
                    aes(x = PC1, y = PC2, label = Variable),
                    size = 4, inherit.aes = FALSE, color = "black") +
    theme_minimal() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = "PCA of Bioclim Variables for Cephalanthera austiniae",
         x = paste0("PC1 (", round(summary(pca_res)$importance[2,1] * 100, 1), "%)"),

         y = paste0("PC2 (", round(summary(pca_res)$importance[2,2] * 100, 1), "%)")) 
