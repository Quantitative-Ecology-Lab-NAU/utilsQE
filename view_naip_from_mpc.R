library(sf)
library(terra)
library(rstac)
library(copc4R)

# draw AOI
aoi <- copc4R::get_aoi()
if (is.null(aoi)) stop("No AOI drawn.")

# keep AOI in lon/lat for STAC search
aoi <- sf::st_transform(aoi, 4326)

# STAC client
stac_obj <- rstac::stac("https://planetarycomputer.microsoft.com/api/stac/v1")

# bbox for STAC search
bb <- sf::st_bbox(aoi)
bbox_vec <- c(
  unname(bb["xmin"]),
  unname(bb["ymin"]),
  unname(bb["xmax"]),
  unname(bb["ymax"])
)

# search NAIP
items <- stac_obj |>
  rstac::stac_search(
    collections = "naip",
    bbox = bbox_vec,
    limit = 100
  ) |>
  rstac::post_request()

items <- rstac::items_sign_planetary_computer(items)

if (length(items$features) == 0) {
  stop("No NAIP imagery found for AOI.")
}

# choose most recent item
item_dates <- vapply(
  items$features,
  function(x) {
    d <- x$properties$datetime
    if (is.null(d)) d <- NA_character_
    d
  },
  character(1)
)

item <- items$features[[order(as.POSIXct(item_dates), decreasing = TRUE)[1]]]
href <- item$assets$image$href

# GDAL settings for remote COG access
Sys.setenv(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_USE_HEAD = "NO"
)

# open remote raster
r <- terra::rast(paste0("/vsicurl/", href))

# project AOI with sf, then convert to terra
aoi_r <- sf::st_transform(aoi, terra::crs(r))
aoi_v <- terra::vect(aoi_r)

# crop and mask by plain extent only
e <- terra::ext(aoi_v)
rc <- terra::crop(r, e, snap = "out")
rc <- terra::mask(rc, e)

# plot
terra::plotRGB(rc, r = 1, g = 2, b = 3, stretch = "lin")
