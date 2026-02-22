###############################################################################
# process_lidar_cloud2trees.R
#
# Complete workflow for processing LiDAR data with topographic feature 
# classification and cloud2trees forest analysis pipeline.
#
# This script:
# 1. Reads LAZ data into a catalog (for one tile for now)
# 2. Classifies ground and normalizes point cloud
# 3. Classifies cliffs and hoodoos (topographic features)
# 4. Runs cloud2trees pipeline for tree detection, segmentation, and metrics
# 5. Outputs tree, forest, canopy fuels, and biomass data
#
# Author: BRCA Analysis
# Date: 2026-02-21
###############################################################################

# =============================================================================
# Load Required Libraries
# =============================================================================

library(lidR)
library(terra)
library(data.table)
library(spanner)
library(RANN)
library(future)
library(sf)
library(ggplot2)

# Install cloud2trees if needed
# remotes::install_github("georgewoolsey/cloud2trees")
library(cloud2trees)

# Source the topographic classification function
source("C:/Users/ajsm/Desktop/BRCA_For_Meag/classify_topo_features.R")

# =============================================================================
# Setup Catalog and Processing Options
# =============================================================================

# Define catalog path - it's one file for now, but chnage this to the folder for the entire catalog
catalog_path <- "//134.114.69.2/QuantEcoData/ALS_Data/Y2019_Bryce_Canyon/2019Lidar/USGS_LPC_UT_Southern_QL1_2018_12SUG9560_LAS_2019.laz"

# Output directory
output_dir <- "C:/Users/ajsm/Desktop/BRCA_For_Meag/outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Intermediate LAS output directories
classified_chunks_dir <- file.path(output_dir, "classified_chunks")
trees_chunks_dir <- file.path(output_dir, "trees_chunks")
dir.create(classified_chunks_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(trees_chunks_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Resolution Parameters
# =============================================================================

dtm_res <- 0.5              # DTM resolution (meters)
chm_res <- 0.5              # Canopy height model resolution (meters)
metrics_res <- 10           # Pixel metrics resolution (meters)

# =============================================================================
# Processing Options
# =============================================================================

write_intermediate_las <- TRUE  # Write intermediate LAS files to disk (reduces memory usage)

# =============================================================================
# Ground Classification Parameters (CSF - tuned for rugged terrain)
# =============================================================================

ground_params <- list(
  sloop_smooth     = FALSE,
  cloth_resolution = 0.25,
  rigidness        = 1L,
  class_threshold  = 0.20,
  iterations       = 800L,
  time_step        = 0.65
)

# =============================================================================
# Topographic Feature Classification Parameters
# =============================================================================

# Parameters for cliff detection
cliff_params <- list(
  slope_thr = 65,              # Slope threshold in degrees
  rescue_h = 0.35,             # Near-ground rescue threshold
  mask_expand_k = 9,           # Focal kernel size
  edge_buffer_m = 8,           # Buffer distance from cliff edges
  edge_max_h_m = 3.0,          # Max height in cliff edge zone
  iso_min_h = 5.0,             # Min height for isolation test
  iso_search_r = 3.0,          # Horizontal search radius
  iso_z_tol = 2.0,             # Vertical tolerance
  iso_min_neighbors = 10,      # Min neighbors required
  max_height = 35              # Max expected canopy height
)

# Parameters for hoodoo detection
hoodoo_params <- list(
  min_height = 6,              # Minimum height for hoodoos
  top_planarity = 0.30,        # Top planarity threshold
  cap_break = 1.20,            # Cap break ratio threshold
  axis_wobble = 0.18,          # Axis wobble threshold
  near_steep_m = 15            # Distance from steep areas
)

# =============================================================================
# cloud2trees Parameters
# =============================================================================

# Tree detection parameters (using local maximum filter)
tree_detection_params <- list(
  ws = function(x) { x * 0.06 + 3 },  # Dynamic window size based on height
  hmin = 2                             # Minimum tree height (meters)
)

# Segmentation parameters
segmentation_params <- list(
  algo = "dalponte2016",         # Segmentation algorithm
  th_tree = 2,                   # Minimum height to be considered tree
  th_seed = 0.45,                # Seed threshold
  th_cr = 0.55,                  # Crown threshold
  max_cr = 10                    # Maximum crown diameter
)

# Metrics to calculate
metrics_params <- list(
  calculate_biomass = TRUE,      # Calculate biomass
  calculate_fuels = TRUE,        # Calculate canopy fuels
  calculate_diversity = TRUE     # Calculate structural diversity
)

# =============================================================================
# Integrated Processing Function
# =============================================================================

#' Complete LiDAR Processing and cloud2trees Pipeline
#' 
#' Processes LAS data through complete workflow:
#' 1. Ground classification
#' 2. Normalization
#' 3. Topographic feature classification (cliffs, hoodoos)
#' 4. Tree detection, segmentation, and metric calculation
#' 
#' @param chunk LAScatalog chunk (passed by catalog_apply)
#' @param ground_params Ground classification parameters
#' @param cliff_params Cliff detection parameters
#' @param hoodoo_params Hoodoo detection parameters
#' @param tree_detection_params Tree detection parameters
#' @param segmentation_params Tree segmentation parameters
#' @param classified_chunks_dir Directory to save classified LAS chunks
#' @param trees_chunks_dir Directory to save tree-segmented LAS chunks
#' @param write_intermediate_las Logical; if TRUE, writes LAS files to disk
#' @return List containing all processed data and metrics
run_cloud2trees_pipeline <- function(chunk, ground_params, cliff_params, hoodoo_params, 
                                      tree_detection_params, segmentation_params,
                                      classified_chunks_dir, trees_chunks_dir,
                                      write_intermediate_las = TRUE) {
  
  # Read the chunk into a las  
  las <- readLAS(chunk)
  
  # Check if chunk has points
  if (is.empty(las)) {
    message("Empty chunk, skipping...")
    return(NULL)
  }
  
  # Get chunk info for logging
  chunk_bbox <- st_bbox(las)
  chunk_center_x <- mean(c(chunk_bbox["xmin"], chunk_bbox["xmax"]))
  chunk_center_y <- mean(c(chunk_bbox["ymin"], chunk_bbox["ymax"]))
  
  cat(
    "\n",
    "================================================================================\n",
    sprintf("=== Processing Chunk: Center (%.1f, %.1f) ===\n", chunk_center_x, chunk_center_y),
    sprintf("=== Extent: X[%.1f - %.1f], Y[%.1f - %.1f] ===\n", 
            chunk_bbox["xmin"], chunk_bbox["xmax"],
            chunk_bbox["ymin"], chunk_bbox["ymax"]),
    "================================================================================\n",
    sprintf("Points: %d\n", npoints(las)),
    sep = ""
  )
  
  # Step 1: Classify ground
  cat("\n1. Classifying ground...\n")
  las <- classify_ground(las, csf(
    sloop_smooth     = ground_params$sloop_smooth,
    cloth_resolution = ground_params$cloth_resolution,
    rigidness        = ground_params$rigidness,
    class_threshold  = ground_params$class_threshold,
    iterations       = ground_params$iterations,
    time_step        = ground_params$time_step
  ))
  
  # Step 2: Build DTM and normalize
  cat("2. Building DTM and normalizing...\n")
  dtm <- rasterize_terrain(las, res = dtm_res, algorithm = tin())
  las_norm <- normalize_height(las, dtm)
  
  # Step 3: Classify topographic features (cliffs and hoodoos)
  cat("3. Classifying topographic features (cliffs, hoodoos)...\n")
  las_classified <- classify_topo_features(
    las_norm,
    features = c("cliff", "hoodoo"),
    classes = list(cliff = 19L, hoodoo = 20L),
    dtm_res = dtm_res,
    cliff_params = cliff_params,
    hoodoo_params = hoodoo_params,
    ncpu = 2,
    verbose = TRUE
  )
  
  # Step 4: Filter out topographic features for vegetation analysis
  cat("4. Filtering vegetation points...\n")
  # Keep only vegetation and ground points (exclude classes 19, 20)
  las_veg <- filter_poi(las_classified, !Classification %in% c(19, 20))
  
  cat("\n=== Running cloud2trees Pipeline ===\n")
  
  # Step 5: Tree detection (Individual Tree Detection - ITD)
  cat("5. Detecting tree tops...\n")
  ttops <- locate_trees(
    las_veg, 
    lmf(ws = tree_detection_params$ws, hmin = tree_detection_params$hmin),
    uniqueness = "bitmerge"
  )
  
  cat("   Detected", nrow(ttops), "trees\n")
  
  # Step 6: Tree segmentation
  cat("6. Segmenting tree crowns...\n")
  las_trees <- segment_trees(
    las_veg,
    dalponte2016(
      chm = rasterize_canopy(las_veg, res = chm_res, algorithm = p2r()),
      treetops = ttops,
      th_tree = segmentation_params$th_tree,
      th_seed = segmentation_params$th_seed,
      th_cr = segmentation_params$th_cr,
      max_cr = segmentation_params$max_cr
    )
  )
  
  # Remove points not assigned to trees
  trees <- filter_poi(las_trees, !is.na(treeID))
  n_trees_segmented <- uniqueN(trees$treeID)
  cat("   Segmented", n_trees_segmented, "trees\n")
  
  # Step 7: Extract tree metrics (points and polygons)
  cat("7. Calculating tree-level metrics...\n")
  tree_metrics <- crown_metrics(trees, .stdtreemetrics)
  tree_metrics <- as.data.table(tree_metrics)
  
  # Add XYZ coordinates
  if (!"X" %in% names(tree_metrics)) {
    tree_metrics[, c("X", "Y", "Z") := .(st_coordinates(geometry)[,1],
                                          st_coordinates(geometry)[,2],
                                          Z)]
  }
  
  # Extract crown polygons (convex hulls)
  cat("   Extracting crown polygons...\n")
  crown_polygons <- crown_metrics(trees, .stdtreemetrics, geom = "convex")
  cat("   Extracted", nrow(crown_polygons), "crown polygons\n")
  
  # Step 8: Calculate biomass (using allometric equations)
  cat("8. Calculating tree biomass...\n")
  # Simple biomass estimation (Jenkins et al. 2003 - generic softwood)
  # Biomass (kg) = exp(β0 + β1 * ln(DBH))
  # DBH estimated from crown area: DBH ≈ 2.6 * sqrt(CrownArea)
  
  if ("convhull_area" %in% names(tree_metrics)) {
    tree_metrics[, DBH_est := 2.6 * sqrt(convhull_area)]
    tree_metrics[, biomass_kg := exp(-2.5356 + 2.4349 * log(pmax(DBH_est, 2.5)))]
    tree_metrics[, carbon_kg := biomass_kg * 0.5]
  }
  
  # Step 9: Calculate canopy metrics and fuels
  cat("9. Calculating canopy and fuel metrics...\n")
  
  # Canopy height model
  # chm <- rasterize_canopy(trees, res = chm_res, algorithm = p2r())
  w <- matrix(1, 3, 3)
  
  fill_na_edge_aware <- function(x, i=5, max_drop=0.35, min_n=3) {
    # If center is valid, keep it (no smoothing of canopy)
    if (!is.na(x[i])) return(x[i])
    
    nb <- x[-i]
    nb <- nb[!is.na(nb)]
    if (length(nb) < min_n) return(NA_real_)
    
    # Use a robust reference (median) and drop neighbors that are much lower
    ref <- median(nb)
    keep <- nb[nb >= (1 - max_drop) * ref]  # e.g., keep within 35% below ref
    if (length(keep) < min_n) keep <- nb    # fallback if too strict
    
    # Robust average: trimmed mean reduces influence of low edge pixels
    mean(keep, trim = 0.2)
  }
  
  chm <- rasterize_canopy(trees, res = chm_res, algorithm = p2r(subcircle = 0.15), pkg = "terra")
  chm[chm < 1.37] <- NA
  
  chm <- terra::focal(chm, w, fun = fill_na_edge_aware, na.policy = "omit")

  # Canopy cover
  canopy_cover <- pixel_metrics(trees, ~mean(Z > 2), res = metrics_res)
  
  # Fuel metrics (height strata)
  fuel_metrics <- pixel_metrics(
    trees,
    ~list(
      surface_fuel = sum(Z >= 0 & Z < 0.5) / .N,
      ladder_fuel = sum(Z >= 0.5 & Z < 4) / .N,
      canopy_fuel = sum(Z >= 4) / .N,
      mean_height = mean(Z),
      max_height = max(Z),
      cbd = sum(Z > 2) / (max(Z) - 2 + 0.001)  # Canopy bulk density proxy
    ),
    res = metrics_res
  )
  
  # Step 10: Forest-level summary statistics
  cat("10. Calculating forest-level statistics...\n")
  forest_summary <- list(
    total_points = npoints(las_classified),
    ground_points = npoints(filter_poi(las_classified, Classification == 2)),
    veg_points = npoints(trees),
    num_trees = n_trees_segmented,
    mean_height = mean(tree_metrics$Z, na.rm = TRUE),
    max_height = max(tree_metrics$Z, na.rm = TRUE),
    mean_crown_area = mean(tree_metrics$convhull_area, na.rm = TRUE),
    total_biomass = sum(tree_metrics$biomass_kg, na.rm = TRUE),
    mean_biomass = mean(tree_metrics$biomass_kg, na.rm = TRUE),
    total_carbon_kg = sum(tree_metrics$carbon_kg, na.rm = TRUE),
    basal_area_m2 = sum(pi * (tree_metrics$DBH_est/200)^2, na.rm = TRUE),
    trees_per_ha = n_trees_segmented / (st_area(st_as_sfc(st_bbox(trees))) / 10000),
    mean_canopy_cover = mean(values(canopy_cover), na.rm = TRUE)
  )
  
  cat(
    "\nForest Summary:\n",
    sprintf("  Trees: %d\n", forest_summary$num_trees),
    sprintf("  Mean Height: %.2f m\n", round(forest_summary$mean_height, 2)),
    sprintf("  Total Biomass: %.0f kg\n", round(forest_summary$total_biomass, 0)),
    sprintf("  Trees/ha: %.1f\n", round(forest_summary$trees_per_ha, 1)),
    sep = ""
  )
  
  # Step 11: Handle intermediate LAS files (write to disk or keep in memory)
  # Get the core (non-buffered) extent from the chunk
  # The 'chunk' object contains the processing extent without buffer
  core_bbox <- st_bbox(chunk)
  
  # Clip outputs to remove buffer (consistent with saved LAS files)
  cat("11. Clipping outputs to remove buffer...\n")
  
  # Clip ttops (sf object) using st_crop
  ttops <- st_crop(ttops, core_bbox)
  
  # Filter tree_metrics based on coordinates (it's a data.table with geometry column)
  # Keep only points within core extent
  tree_metrics <- tree_metrics[
    X >= core_bbox["xmin"] & X <= core_bbox["xmax"] &
    Y >= core_bbox["ymin"] & Y <= core_bbox["ymax"]
  ]
  
  # Filter crown polygons to match treeIDs from filtered tree_metrics
  # Only keep polygons whose treeID exists in the core extent tree_metrics
  crown_polygons <- crown_polygons[crown_polygons$treeID %in% tree_metrics$treeID, ]
  
  # Crop rasters to core extent using terra::crop with terra::ext
  dtm <- terra::crop(dtm, terra::ext(core_bbox["xmin"], core_bbox["xmax"], 
                                      core_bbox["ymin"], core_bbox["ymax"]))
  chm <- terra::crop(chm, terra::ext(core_bbox["xmin"], core_bbox["xmax"], 
                                      core_bbox["ymin"], core_bbox["ymax"]))
  canopy_cover <- terra::crop(canopy_cover, terra::ext(core_bbox["xmin"], core_bbox["xmax"], 
                                                        core_bbox["ymin"], core_bbox["ymax"]))
  fuel_metrics <- terra::crop(fuel_metrics, terra::ext(core_bbox["xmin"], core_bbox["xmax"], 
                                                        core_bbox["ymin"], core_bbox["ymax"]))
  
  cat("   Clipped SF objects and rasters to core extent\n")
  
  las_classified_output <- NULL
  las_trees_output <- NULL
  
  if (write_intermediate_las) {
    cat("12. Writing intermediate LAS files to disk...\n")
    # Create unique filenames based on chunk center coordinates
    chunk_id <- sprintf("%.0f_%.0f", chunk_center_x, chunk_center_y)
    
    # Clip LAS objects to remove buffer before saving
    cat("   Clipping LAS files to core extent...\n")
    las_classified_core <- clip_rectangle(las_classified, 
                                          xleft = core_bbox["xmin"], 
                                          ybottom = core_bbox["ymin"],
                                          xright = core_bbox["xmax"], 
                                          ytop = core_bbox["ymax"])
    
    las_trees_core <- clip_rectangle(las_trees, 
                                     xleft = core_bbox["xmin"], 
                                     ybottom = core_bbox["ymin"],
                                     xright = core_bbox["xmax"], 
                                     ytop = core_bbox["ymax"])
    
    # Write classified LAS (compressed as LAZ, without buffer)
    las_classified_output <- file.path(classified_chunks_dir, paste0("classified_", chunk_id, ".laz"))
    writeLAS(las_classified_core, las_classified_output)
    cat("   Saved classified LAS:", basename(las_classified_output), "\n")
    
    # Write segmented trees LAS (compressed as LAZ, without buffer)
    las_trees_output <- file.path(trees_chunks_dir, paste0("trees_", chunk_id, ".laz"))
    writeLAS(las_trees_core, las_trees_output)
    cat("   Saved trees LAS:", basename(las_trees_output), "\n")
  } else {
    cat("12. Keeping LAS objects in memory (write_intermediate_las = FALSE)...\n")
    las_classified_output <- las_classified
    las_trees_output <- las_trees
  }
  
  # Return comprehensive results
  # Wrap SpatRaster objects for safe serialization across parallel workers
  return(list(
    las_classified_file = if(write_intermediate_las) las_classified_output else NULL,
    las_trees_file = if(write_intermediate_las) las_trees_output else NULL,
    las_classified = if(!write_intermediate_las) las_classified_output else NULL,
    las_trees = if(!write_intermediate_las) las_trees_output else NULL,
    tree_tops = ttops,
    tree_metrics = tree_metrics,
    crown_polygons = crown_polygons,
    forest_summary = forest_summary,
    dtm = terra::wrap(dtm),
    chm = terra::wrap(chm),
    canopy_cover = terra::wrap(canopy_cover),
    fuel_metrics = terra::wrap(fuel_metrics)
  ))
}

# =============================================================================
# Main Processing Workflow
# =============================================================================

cat("\n",
  "================================================================================\n",
  "  LiDAR Processing Workflow with Topographic Classification and cloud2trees  \n",
  "================================================================================\n",
  "\n", sep = "")

# Create catalog
cat("\nReading catalog...\n")
ctg <- readLAScatalog(catalog_path)

cat("\nCatalog summary:\n")
print(ctg)

# Configure catalog processing options
opt_chunk_buffer(ctg) <- 20        # 20m buffer for edge effects
opt_chunk_size(ctg) <- 300         # Process 300m sq chunks (0 = no chunking)
opt_output_files(ctg) <- ""        # In-memory processing
opt_progress(ctg) <- TRUE          # Show progress bar
opt_stop_early(ctg) <- TRUE       # Process all chunks

# Configure catalog filter to drop noise and withheld points
opt_filter(ctg) <- "-drop_class 7 18"

cat("\nCatalog processing options configured:\n",
  sprintf("  Chunk size: %d x %d m\n", opt_chunk_size(ctg), opt_chunk_size(ctg)),
  sprintf("  Chunk buffer: %d m\n", opt_chunk_buffer(ctg)),
  sprintf("  Filter: %s\n", opt_filter(ctg)),
  sprintf("  Progress: %s (will display chunk plot)\n", opt_progress(ctg)),
  sep = "")

# Set up parallel processing with future
plan(multisession, workers = 6)    # Adjust workers based on your CPU cores
set_lidr_threads(1)                # Threads per worker
cat("\nSetting up parallel processing...\n",
    "  Workers: 6\n",
    "  Threads per worker: 1\n",
    sep = "")

# Process using catalog_apply
cat("\n=== Starting Catalog Processing ===\n", 
    "Processing will show chunk progress below...\n\n",
    sep = "")

results <- catalog_apply(ctg, run_cloud2trees_pipeline, 
                         ground_params = ground_params,
                         cliff_params = cliff_params,
                         hoodoo_params = hoodoo_params,
                         tree_detection_params = tree_detection_params,
                         segmentation_params = segmentation_params,
                         classified_chunks_dir = classified_chunks_dir,
                         trees_chunks_dir = trees_chunks_dir,
                         write_intermediate_las = write_intermediate_las)

cat("\n=== Catalog Processing Complete ===\n")

# =============================================================================
# Merge Results from Chunks
# =============================================================================

cat("\n=== Merging results from chunks ===\n")

# Note: LAS objects are NOT merged to save memory.
# Individual chunks are saved in classified_chunks/ and trees_chunks/ folders.

# Merge tree tops (SF objects)
cat("  Merging tree tops...\n")
tree_tops_list <- lapply(results, function(x) x$tree_tops)
tree_tops_list <- tree_tops_list[!sapply(tree_tops_list, is.null)]
if (length(tree_tops_list) > 0) {
  tree_tops_merged <- do.call(rbind, tree_tops_list)
} else {
  stop("No tree tops detected in any chunks")
}

# Merge tree metrics (data.table)
cat("  Merging tree metrics...\n")
tree_metrics_list <- lapply(results, function(x) x$tree_metrics)
tree_metrics_list <- tree_metrics_list[!sapply(tree_metrics_list, is.null)]
tree_metrics_list <- tree_metrics_list[sapply(tree_metrics_list, nrow) > 0]
if (length(tree_metrics_list) > 0) {
  tree_metrics_merged <- rbindlist(tree_metrics_list, fill = TRUE)
} else {
  stop("No tree metrics found in any chunks")
}

# Merge crown polygons (SF objects)
cat("  Merging crown polygons...\n")
crown_polygons_list <- lapply(results, function(x) x$crown_polygons)
crown_polygons_list <- crown_polygons_list[!sapply(crown_polygons_list, is.null)]
if (length(crown_polygons_list) > 0) {
  crown_polygons_merged <- do.call(rbind, crown_polygons_list)
  cat(sprintf("    Merged %d crown polygons\n", nrow(crown_polygons_merged)))
} else {
  stop("No crown polygons extracted from any chunks")
}

# Merge rasters using terra::merge
cat("  Merging rasters...\n")
# Unwrap rasters first (they were wrapped for serialization)
dtm_list <- lapply(results, function(x) terra::rast(x$dtm))
dtm_list <- dtm_list[!sapply(dtm_list, is.null)]
if (length(dtm_list) == 1) {
  dtm_merged <- dtm_list[[1]]
  cat("    DTM: single raster, no merge needed\n")
} else if (length(dtm_list) > 1) {
  cat(sprintf("    DTM: merging %d rasters...\n", length(dtm_list)))
  dtm_merged <- tryCatch({
    do.call(terra::merge, dtm_list)
  }, error = function(e) {
    cat("    Warning: merge failed, trying mosaic...\n")
    do.call(terra::mosaic, dtm_list)
  })
} else {
  stop("No DTM rasters found")
}

chm_list <- lapply(results, function(x) terra::rast(x$chm))
chm_list <- chm_list[!sapply(chm_list, is.null)]
if (length(chm_list) == 1) {
  chm_merged <- chm_list[[1]]
  cat("    CHM: single raster, no merge needed\n")
} else if (length(chm_list) > 1) {
  cat(sprintf("    CHM: merging %d rasters...\n", length(chm_list)))
  chm_merged <- tryCatch({
    do.call(terra::merge, chm_list)
  }, error = function(e) {
    cat("    Warning: merge failed, trying mosaic...\n")
    do.call(terra::mosaic, chm_list)
  })
} else {
  stop("No CHM rasters found")
}

canopy_cover_list <- lapply(results, function(x) terra::rast(x$canopy_cover))
canopy_cover_list <- canopy_cover_list[!sapply(canopy_cover_list, is.null)]
if (length(canopy_cover_list) == 1) {
  canopy_cover_merged <- canopy_cover_list[[1]]
  cat("    Canopy cover: single raster, no merge needed\n")
} else if (length(canopy_cover_list) > 1) {
  cat(sprintf("    Canopy cover: merging %d rasters...\n", length(canopy_cover_list)))
  canopy_cover_merged <- tryCatch({
    do.call(terra::merge, canopy_cover_list)
  }, error = function(e) {
    cat("    Warning: merge failed, trying mosaic...\n")
    do.call(terra::mosaic, canopy_cover_list)
  })
} else {
  stop("No canopy cover rasters found")
}

fuel_metrics_list <- lapply(results, function(x) terra::rast(x$fuel_metrics))
fuel_metrics_list <- fuel_metrics_list[!sapply(fuel_metrics_list, is.null)]
if (length(fuel_metrics_list) == 1) {
  fuel_metrics_merged <- fuel_metrics_list[[1]]
  cat("    Fuel metrics: single raster, no merge needed\n")
} else if (length(fuel_metrics_list) > 1) {
  cat(sprintf("    Fuel metrics: merging %d rasters...\n", length(fuel_metrics_list)))
  fuel_metrics_merged <- tryCatch({
    do.call(terra::merge, fuel_metrics_list)
  }, error = function(e) {
    cat("    Warning: merge failed, trying mosaic...\n")
    do.call(terra::mosaic, fuel_metrics_list)
  })
} else {
  stop("No fuel metrics rasters found")
}

# Aggregate forest summary stats
cat("  Aggregating forest summary...\n")
forest_summaries <- lapply(results, function(x) x$forest_summary)
forest_summaries <- forest_summaries[!sapply(forest_summaries, is.null)]

# Safe extraction function that handles NULL values
safe_sum <- function(summaries, field) {
  vals <- sapply(summaries, function(x) if(is.null(x[[field]])) 0 else x[[field]])
  sum(as.numeric(vals), na.rm = TRUE)
}
safe_mean <- function(summaries, field) {
  vals <- sapply(summaries, function(x) if(is.null(x[[field]])) NA else x[[field]])
  mean(as.numeric(vals), na.rm = TRUE)
}
safe_max <- function(summaries, field) {
  vals <- sapply(summaries, function(x) if(is.null(x[[field]])) NA else x[[field]])
  max(as.numeric(vals), na.rm = TRUE)
}

forest_summary_merged <- list(
  total_points = safe_sum(forest_summaries, "total_points"),
  ground_points = safe_sum(forest_summaries, "ground_points"),
  veg_points = safe_sum(forest_summaries, "veg_points"),
  num_trees = safe_sum(forest_summaries, "num_trees"),
  mean_height = safe_mean(forest_summaries, "mean_height"),
  max_height = safe_max(forest_summaries, "max_height"),
  mean_canopy_cover = safe_mean(forest_summaries, "mean_canopy_cover"),
  mean_biomass = safe_mean(forest_summaries, "mean_biomass"),
  total_biomass = safe_sum(forest_summaries, "total_biomass")
)

cat("  Merge complete!\n")

# =============================================================================
# Save Outputs
# =============================================================================

cat("\n=== Saving outputs ===\n", sep = "")

# Note: Final merged LAS files are NOT saved to reduce disk usage.
# Intermediate chunks are available in classified_chunks/ and trees_chunks/ folders.

# Save tree metrics as CSV
output_tree_metrics <- file.path(output_dir, "brca_tree_metrics.csv")
fwrite(tree_metrics_merged, output_tree_metrics)
cat("  Saved tree metrics CSV:", output_tree_metrics, "\n")

# Save spatial data to single GeoPackage with multiple layers
output_gpkg <- file.path(output_dir, "brca_spatial_outputs.gpkg")

# Delete existing gpkg if it exists
if (file.exists(output_gpkg)) {
  file.remove(output_gpkg)
}

# Layer 1: Tree tops (3D points)
st_write(tree_tops_merged, output_gpkg, layer = "tree_tops", delete_dsn = FALSE, quiet = TRUE)
cat("  Saved to GeoPackage - Layer 'tree_tops':", output_gpkg, "\n")

# Layer 2: Tree metrics as points (if it has X,Y coordinates)
if ("X" %in% colnames(tree_metrics_merged) && "Y" %in% colnames(tree_metrics_merged)) {
  tree_metrics_sf <- st_as_sf(tree_metrics_merged, 
                               coords = c("X", "Y"), 
                               crs = st_crs(tree_tops_merged),
                               remove = FALSE)
  st_write(tree_metrics_sf, output_gpkg, layer = "tree_metrics", append = TRUE, quiet = TRUE)
  cat("  Saved to GeoPackage - Layer 'tree_metrics':", output_gpkg, "\n")
}

# Layer 3: Crown polygons
st_write(crown_polygons_merged, output_gpkg, layer = "crown_polygons", append = TRUE, quiet = TRUE)
cat("  Saved to GeoPackage - Layer 'crown_polygons':", output_gpkg, "\n")

cat("  GeoPackage complete:", output_gpkg, "\n")

# Save rasters
output_dtm <- file.path(output_dir, "brca_dtm.tif")
writeRaster(dtm_merged, output_dtm, overwrite = TRUE)
cat("  Saved DTM:", output_dtm, "\n")

output_chm <- file.path(output_dir, "brca_chm.tif")
writeRaster(chm_merged, output_chm, overwrite = TRUE)
cat("  Saved CHM:", output_chm, "\n")

output_canopy_cover <- file.path(output_dir, "brca_canopy_cover.tif")
writeRaster(canopy_cover_merged, output_canopy_cover, overwrite = TRUE)
cat("  Saved canopy cover:", output_canopy_cover, "\n")

output_fuel_metrics <- file.path(output_dir, "brca_fuel_metrics.tif")
writeRaster(fuel_metrics_merged, output_fuel_metrics, overwrite = TRUE)
cat("  Saved fuel metrics:", output_fuel_metrics, "\n")

# Save forest summary as JSON
output_forest_summary <- file.path(output_dir, "brca_forest_summary.json")
jsonlite::write_json(forest_summary_merged, output_forest_summary, 
                     pretty = TRUE, auto_unbox = TRUE)
cat("  Saved forest summary:", output_forest_summary, "\n")

# =============================================================================
# Visualizations
# =============================================================================

cat("\n=== Creating visualizations ===\n", sep = "")

# 1. Canopy Height Model (CHM)
cat("  Creating CHM plot...\n")
chm_df <- as.data.frame(chm_merged, xy = TRUE)
colnames(chm_df) <- c("x", "y", "height")

p1 <- ggplot(chm_df, aes(x = x, y = y, fill = height)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", name = "Height (m)", na.value = "transparent") +
  coord_equal() +
  labs(title = "Canopy Height Model",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "plot_chm.png"), p1, width = 10, height = 8, dpi = 300)
cat("  Saved CHM plot\n")

# 2. Canopy Cover
cat("  Creating canopy cover plot...\n")
canopy_df <- as.data.frame(canopy_cover_merged, xy = TRUE)
colnames(canopy_df) <- c("x", "y", "cover")

p2 <- ggplot(canopy_df, aes(x = x, y = y, fill = cover)) +
  geom_raster() +
  scale_fill_viridis_c(option = "mako", name = "Cover (%)", na.value = "transparent") +
  coord_equal() +
  labs(title = "Canopy Cover",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "plot_canopy_cover.png"), p2, width = 10, height = 8, dpi = 300)
cat("  Saved canopy cover plot\n")

# 3. Fuel Metrics
cat("  Creating fuel metrics plot...\n")
fuel_df <- as.data.frame(fuel_metrics_merged, xy = TRUE)
colnames(fuel_df) <- c("x", "y", "fuel")

p4 <- ggplot(fuel_df, aes(x = x, y = y, fill = fuel)) +
  geom_raster() +
  scale_fill_viridis_c(option = "inferno", name = "Fuel Load", na.value = "transparent") +
  coord_equal() +
  labs(title = "Fuel Metrics",
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "plot_fuel_metrics.png"), p4, width = 10, height = 8, dpi = 300)
cat("  Saved fuel metrics plot\n")

# 4. Tree Tops
cat("  Creating tree tops plot...\n")
tree_tops_coords <- st_coordinates(tree_tops_merged)
tree_tops_attrs <- st_drop_geometry(tree_tops_merged)
# Keep Z from attributes if it exists, otherwise use Z from coordinates
if (!"Z" %in% colnames(tree_tops_attrs) && "Z" %in% colnames(tree_tops_coords)) {
  # Z is in coordinates (3D points), add it to data frame
  tree_tops_data <- cbind(as.data.frame(tree_tops_coords), 
                          tree_tops_attrs[, !colnames(tree_tops_attrs) %in% c("X", "Y"), drop = FALSE])
} else {
  # Z is in attributes, remove X, Y from attributes to avoid duplicates
  tree_tops_data <- cbind(as.data.frame(tree_tops_coords[, c("X", "Y"), drop = FALSE]), 
                          tree_tops_attrs[, !colnames(tree_tops_attrs) %in% c("X", "Y"), drop = FALSE])
}

p5 <- ggplot(tree_tops_data, aes(x = X, y = Y, color = Z)) +
  geom_point(size = 1.5, alpha = 0.6) +
  scale_color_viridis_c(option = "turbo", name = "Height (m)") +
  coord_equal() +
  labs(title = "Detected Tree Tops",
       subtitle = paste0("Total trees: ", nrow(tree_tops_merged)),
       x = "Easting (m)", y = "Northing (m)") +
  theme_minimal() +
  theme(legend.position = "right")

ggsave(file.path(output_dir, "plot_tree_tops.png"), p5, width = 10, height = 8, dpi = 300)
cat("  Saved tree tops plot\n")

# 5. Tree Metrics (colored by biomass or height)
cat("  Creating tree metrics plot...\n")
if ("biomass_kg" %in% colnames(tree_metrics_merged)) {
  p6 <- ggplot(tree_metrics_merged, aes(x = X, y = Y, color = biomass_kg, size = Z)) +
    geom_point(alpha = 0.6) +
    scale_color_viridis_c(option = "rocket", name = "Biomass (kg)", trans = "log10") +
    scale_size_continuous(name = "Height (m)", range = c(0.5, 3)) +
    coord_equal() +
    labs(title = "Tree Metrics - Biomass",
         subtitle = paste0("N = ", nrow(tree_metrics_merged), " trees"),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "plot_tree_biomass.png"), p6, width = 10, height = 8, dpi = 300)
  cat("  Saved tree biomass plot\n")
} else {
  p6 <- ggplot(tree_metrics_merged, aes(x = X, y = Y, color = Z)) +
    geom_point(size = 2, alpha = 0.6) +
    scale_color_viridis_c(option = "rocket", name = "Height (m)") +
    coord_equal() +
    labs(title = "Tree Metrics - Height",
         subtitle = paste0("N = ", nrow(tree_metrics_merged), " trees"),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(legend.position = "right")
  
  ggsave(file.path(output_dir, "plot_tree_height.png"), p6, width = 10, height = 8, dpi = 300)
  cat("  Saved tree height plot\n")
}

# 6. Tree Crown Polygons (colored by hull area)
cat("  Creating tree crown polygons plot...\n")
# Extract the hull area from crown_polygons_merged
if ("convhull_area" %in% colnames(crown_polygons_merged)) {
  p7 <- ggplot(crown_polygons_merged) +
    geom_sf(aes(fill = convhull_area), color = "gray40", size = 0.1, alpha = 0.7) +
    scale_fill_viridis_c(option = "plasma", name = "Crown Area (m²)", 
                         trans = "log10",
                         breaks = c(1, 5, 10, 25, 50, 100)) +
    coord_sf() +
    labs(title = "Tree Crown Polygons",
         subtitle = paste0("N = ", nrow(crown_polygons_merged), " trees"),
         x = "Easting (m)", y = "Northing (m)") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid = element_line(color = "gray90"))
  
  ggsave(file.path(output_dir, "plot_crown_polygons.png"), p7, width = 10, height = 8, dpi = 300)
  cat("  Saved tree crown polygons plot\n")
} else {
  cat("  Warning: convhull_area not found in crown_polygons, skipping crown polygon plot\n")
}

cat("  All visualizations saved to:", output_dir, "\n")

# =============================================================================
# Complete!
# =============================================================================

# Clean up
plan(sequential)  # Reset parallel processing

