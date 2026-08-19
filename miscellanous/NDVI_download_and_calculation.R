#install.packages(c("rsi", "sf", "terra"))

library(rsi)
library(sf)
library(terra)
#### setwd
setwd("/Users/gabri/Library/CloudStorage/OneDrive-UniversityofArizona/mice/Suzuki_meeting_presentation/")

#### read files in metadata
metadata = readRDS("metadata_clean.RDS")

metadata$x_coord
metadata$y_coord

# 1) Put your coordinates into a data frame
pts <- data.frame(
  x = metadata$x_coord,   # longitude
  y = metadata$y_coord    # latitude
)

# 2) Convert to sf points in WGS84
pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = 4326)

# 3) Make a bounding box around all points
aoi <- st_as_sfc(st_bbox(pts_sf))

aoi


#pt <- st_sfc(
#  st_point(c(metadata$x_coord[1], metadata$y_coord[1])),
#  crs = 4326
#)

# Buffer in a projected CRS chosen from the point
# Create AOI in a projected CRS
#aoi_proj <- st_buffer(st_transform(pt, 26912), 1000)

# Keep it projected when calling the download function
#landsat_file <- get_landsat_imagery(
#  aoi = aoi_proj,
#  start_date = "2024-06-01",
#  end_date   = "2024-06-30",
#  mask_band = "qa_pixel",
#  output_filename = "landsat_30m.tif"
#)

#r <- rast(landsat_file)
#plot(r[[1]])


#r
#nlyr(r)
#dim(r)
#ext(r)
#crs(r)

#plot(r[[1]])
#plotRGB(r, r = "R", g = "G", b = "B", stretch = "lin")

#file.info(landsat_file)$size / 1024^2




pts_sf <- st_as_sf(
  metadata,
  coords = c("x_coord", "y_coord"),
  crs = 4326,
  remove = FALSE
)

aoi_ll <- st_as_sfc(st_bbox(pts_sf))

aoi_proj <- st_transform(aoi_ll, 26912)
aoi_proj <- st_buffer(aoi_proj, 1000)


### 2025 ------------------
# Keep it projected when calling the download function 
landsat_file <- get_landsat_imagery(
  aoi = aoi_proj,
  start_date = "2025-06-01",
  end_date   = "2025-06-30",
  mask_band = "qa_pixel",
  output_filename = "landsat_30m.tif"
)

r <- rast(landsat_file)
plot(r[[1]])


plot(r[[1]])
plotRGB(r, r = "R", g = "G", b = "B", stretch = "lin")

file.info(landsat_file)$size / 1024^2


ndvi <- (r[["N"]] - r[["R"]]) / (r[["N"]] + r[["R"]])


ndvi_dir <- "~/landsat_ndvi"
dir.create(ndvi_dir, showWarnings = FALSE, recursive = TRUE)

writeRaster(
  ndvi,
  file.path(ndvi_dir, "ndvi_2025_06.tif"),
  overwrite = TRUE
)






### 2024 ------------------
# Keep it projected when calling the download function 
landsat_file <- get_landsat_imagery(
  aoi = aoi_proj,
  start_date = "2024-06-01",
  end_date   = "2024-06-30",
  mask_band = "qa_pixel",
  output_filename = "landsat_30m.tif"
)

r <- rast(landsat_file)
plot(r[[1]])
plotRGB(r, r = "R", g = "G", b = "B", stretch = "lin")

file.info(landsat_file)$size / 1024^2


ndvi <- (r[["N"]] - r[["R"]]) / (r[["N"]] + r[["R"]])


ndvi_dir <- "~/landsat_ndvi"
dir.create(ndvi_dir, showWarnings = FALSE, recursive = TRUE)

writeRaster(
  ndvi,
  file.path(ndvi_dir, "ndvi_2024_06.tif"),
  overwrite = TRUE
)



### 2023 ------------------
# Keep it projected when calling the download function 
landsat_file <- get_landsat_imagery(
  aoi = aoi_proj,
  start_date = "2023-06-01",
  end_date   = "2023-06-30",
  mask_band = "qa_pixel",
  output_filename = "landsat_30m.tif"
)

r <- rast(landsat_file)
plot(r[[1]])
plotRGB(r, r = "R", g = "G", b = "B", stretch = "lin")

file.info(landsat_file)$size / 1024^2


ndvi <- (r[["N"]] - r[["R"]]) / (r[["N"]] + r[["R"]])


ndvi_dir <- "~/landsat_ndvi"

writeRaster(
  ndvi,
  file.path(ndvi_dir, "ndvi_2023_06.tif"),
  overwrite = TRUE
)




### 2022 ------------------
# Keep it projected when calling the download function 
landsat_file <- get_landsat_imagery(
  aoi = aoi_proj,
  start_date = "2022-06-01",
  end_date   = "2022-06-30",
  mask_band = "qa_pixel",
  output_filename = "landsat_30m.tif"
)

r <- rast(landsat_file)
plot(r[[1]])
plotRGB(r, r = "R", g = "G", b = "B", stretch = "lin")

file.info(landsat_file)$size / 1024^2


ndvi <- (r[["N"]] - r[["R"]]) / (r[["N"]] + r[["R"]])


ndvi_dir <- "~/landsat_ndvi"

writeRaster(
  ndvi,
  file.path(ndvi_dir, "ndvi_2022_06.tif"),
  overwrite = TRUE
)





ndvi_dir <- "~/landsat_ndvi"

# Read the four NDVI rasters
ndvi_files <- file.path(
  ndvi_dir,
  c("ndvi_2022_06.tif", "ndvi_2023_06.tif", "ndvi_2024_06.tif", "ndvi_2025_06.tif")
)

template <- rast(ndvi_files[1])

ndvi_aligned <- lapply(ndvi_files, function(f) {
  r <- rast(f)
  
  # force all rasters onto the template grid
  r <- resample(r, template, method = "bilinear")
  
  r
})

ndvi_stack <- rast(ndvi_aligned)

ndvi_mean <- app(ndvi_stack, mean, na.rm = TRUE)

writeRaster(
  ndvi_mean,
  file.path(ndvi_dir, "ndvi_mean_2022_2025_june.tif"),
  overwrite = TRUE
)

plot(ndvi_mean)

library(terra)

# Convert sf points to terra vector
pts_vect <- vect(pts_sf)

# Extract mean NDVI at each point
ndvi_vals <- extract(ndvi_mean, pts_vect)

head(ndvi_vals)
metadata$mean_ndvi_2022_2025 <- ndvi_vals[,2]

plot(metadata$elevation ~ metadata$mean_ndvi_2022_2025)

saveRDS(metadata , "metadata_clean_1.RDS")

