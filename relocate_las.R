#' Relocate a LAS object to a new UTM coordinate
#'
#' This function takes a LAS object or file path, translates it to a new location 
#' specified by target UTM coordinates, properly updates the header offsets, 
#' quantizes the coordinates, and returns a valid LAS object.
#'
#' @param las Either a LAS object or a character string path to a LAS/LAZ file
#' @param target_easting Numeric. Target UTM easting coordinate
#' @param target_northing Numeric. Target UTM northing coordinate
#' @param target_epsg Integer. EPSG code for the target CRS (default: NULL, uses original CRS)
#' @param reference_point Character. Which point to use as reference: "center" (default),
#'   "min", "max", or "bbox_center" for the bounding box center
#' @param select Character. Filter string for selective reading if path provided 
#'   (default: "*" for all points)
#' @param filter Character. Additional filter string for reading (default: "")
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A LAS object relocated to the new coordinates with updated header
#'
#' @details
#' This function performs the following operations:
#' 1. Reads the input LAS file (if path provided)
#' 2. Computes reference point based on current coordinates
#' 3. Calculates translation deltas
#' 4. Updates header offsets to be near the target domain (best practice)
#' 5. Translates point coordinates
#' 6. Sets/validates CRS if EPSG is provided
#' 7. Quantizes coordinates to ensure LAS spec compliance
#' 8. Updates header with new bounding box and metadata
#'
#' @examples
#' \dontrun{
#' # From file path
#' las_relocated <- relocate_las(
#'   las = "input.las",
#'   target_easting = 394026,
#'   target_northing = 4162359,
#'   target_epsg = 26912,
#'   reference_point = "center"
#' )
#' 
#' # Write to file
#' writeLAS(las_relocated, "output.laz")
#' 
#' # From LAS object
#' las <- readLAS("input.las")
#' las_relocated <- relocate_las(las, 394026, 4162359, target_epsg = 26912)
#' writeLAS(las_relocated, "output.laz")
#' }
#'
#' @export
relocate_las <- function(las,
                         target_easting,
                         target_northing,
                         target_epsg = NULL,
                         reference_point = c("center", "min", "max", "bbox_center"),
                         select = "*",
                         filter = "",
                         verbose = TRUE) {
  
  # Load required packages
  if (!requireNamespace("lidR", quietly = TRUE)) {
    stop("Package 'lidR' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required but not installed.", call. = FALSE)
  }
  
  # Validate inputs
  reference_point <- match.arg(reference_point)
  
  if (!is.numeric(target_easting) || !is.numeric(target_northing)) {
    stop("target_easting and target_northing must be numeric", call. = FALSE)
  }
  
  # Read LAS file if path is provided
  if (is.character(las)) {
    if (!file.exists(las)) {
      stop("Input file does not exist: ", las, call. = FALSE)
    }
    if (verbose) message("Reading LAS file: ", las)
    las <- lidR::readLAS(las, select = select, filter = filter)
  } else if (!inherits(las, "LAS")) {
    stop("'las' must be either a file path or a LAS object", call. = FALSE)
  }
  
  # Check for empty LAS
  if (lidR::is.empty(las)) {
    stop("LAS is empty or could not be read", call. = FALSE)
  }
  
  if (verbose) {
    message("  Points: ", nrow(las@data))
    message("  Original extent: X[", round(min(las$X), 2), ", ", round(max(las$X), 2), 
            "], Y[", round(min(las$Y), 2), ", ", round(max(las$Y), 2), "]")
  }
  
  # Calculate reference point based on current coordinates
  x_ref <- switch(reference_point,
                  center = mean(las@data$X),
                  min = min(las@data$X),
                  max = max(las@data$X),
                  bbox_center = (min(las@data$X) + max(las@data$X)) / 2)
  
  y_ref <- switch(reference_point,
                  center = mean(las@data$Y),
                  min = min(las@data$Y),
                  max = max(las@data$Y),
                  bbox_center = (min(las@data$Y) + max(las@data$Y)) / 2)
  
  if (verbose) {
    message("  Reference point (", reference_point, "): X=", round(x_ref, 2), 
            ", Y=", round(y_ref, 2))
  }
  
  # Calculate translation deltas
  dx <- target_easting - x_ref
  dy <- target_northing - y_ref
  
  if (verbose) {
    message("Translation: dX=", round(dx, 2), ", dY=", round(dy, 2))
  }
  
  # Update header offsets to be near the target domain (best practice)
  # Using floor of target coordinates as offsets
  las@header@PHB[["X offset"]] <- floor(target_easting)
  las@header@PHB[["Y offset"]] <- floor(target_northing)
  
  if (verbose) {
    message("Updated header offsets: X offset=", las@header@PHB[["X offset"]], 
            ", Y offset=", las@header@PHB[["Y offset"]])
  }
  
  # Translate point coordinates
  las@data$X <- las@data$X + dx
  las@data$Y <- las@data$Y + dy
  
  # Set or validate CRS
  if (!is.null(target_epsg)) {
    if (verbose) message("Setting CRS to EPSG:", target_epsg)
    sf::st_crs(las) <- sf::st_crs(target_epsg)
  } else {
    if (is.na(sf::st_crs(las))) {
      warning("No CRS set and no target_epsg provided. Output LAS will have no CRS.", 
              call. = FALSE)
    } else {
      if (verbose) message("Keeping original CRS: ", sf::st_crs(las)$input)
    }
  }
  
  # Quantize coordinates to ensure LAS specification compliance
  # This updates the coordinates to match scale/offset precision requirements
  if (verbose) message("Quantizing coordinates...")
  las <- lidR::las_quantize(las, by_reference = FALSE)
  
  # Update header with new bounding box and point statistics
  if (verbose) message("Updating header metadata...")
  las <- lidR::las_update(las)
  
  if (verbose) {
    message("  New extent: X[", round(min(las$X), 2), ", ", round(max(las$X), 2), 
            "], Y[", round(min(las$Y), 2), ", ", round(max(las$Y), 2), "]")
    message("Successfully relocated LAS object!")
  }
  
  return(las)
}
