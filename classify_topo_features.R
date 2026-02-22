###############################################################################
# classify_topo_features.R
#
# Standalone function to identify and classify topographic features (cliffs, 
# hoodoos, depressions) in a normalized LAS object.
#
# Dependencies: lidR, terra, data.table, spanner, RANN
###############################################################################

#' Classify Topographic Features in Normalized LAS
#'
#' Identifies topographic features (cliffs, hoodoos, depressions) in a 
#' normalized LAS object and classifies them into unused ASPRS classes 
#' (19-21) for masking from vegetation analysis.
#'
#' @param las_norm A normalized LAS object (must have Z normalized and Zref column)
#' @param features Character vector of features to classify. Options: "cliff", "hoodoo", "depression"
#' @param classes Named list of ASPRS classes for each feature (default: list(cliff = 19, hoodoo = 20, depression = 21))
#' @param dtm_res DTM resolution in meters (default: 0.5)
#' @param cliff_params List of cliff detection parameters (see details)
#' @param hoodoo_params List of hoodoo detection parameters (see details)
#' @param depression_params List of depression detection parameters (see details)
#' @param ncpu Number of CPUs for parallel processing (default: 8)
#' @param verbose Print progress messages (default: TRUE)
#' 
#' @details 
#' Parameter lists for each feature type:
#' 
#' **cliff_params** (all optional):
#' \itemize{
#'   \item slope_thr: Slope threshold in degrees (default: 65)
#'   \item rescue_h: Near-ground rescue threshold in meters (default: 0.35)
#'   \item mask_expand_k: Focal kernel size for expanding steep mask (default: 9)
#'   \item edge_buffer_m: Buffer distance from cliff edges in meters (default: 8)
#'   \item edge_max_h_m: Max plausible height in cliff edge zone (default: 3.0)
#'   \item iso_min_h: Minimum height to test for isolated points (default: 5.0)
#'   \item iso_search_r: Horizontal search radius for isolation test (default: 3.0)
#'   \item iso_z_tol: Vertical tolerance for neighbor similarity (default: 2.0)
#'   \item iso_min_neighbors: Minimum neighbors required (default: 10)
#'   \item max_height: Maximum expected canopy height (default: 35)
#' }
#' 
#' **hoodoo_params** (all optional):
#' \itemize{
#'   \item min_height: Minimum height for hoodoo detection (default: 6)
#'   \item top_planarity: Top planarity threshold (default: 0.30)
#'   \item cap_break: Cap break ratio threshold (default: 1.20)
#'   \item axis_wobble: Axis wobble threshold (default: 0.18)
#'   \item near_steep_m: Distance from steep areas to search (default: 15)
#' }
#' 
#' **depression_params** (all optional):
#' \itemize{
#'   \item min_depth: Minimum depression depth in meters (default: 0.3)
#'   \item max_depth: Maximum depression depth in meters (default: 5.0)
#'   \item min_area: Minimum depression area in m² (default: 2.0)
#'   \item expand_cells: Number of cells to expand depression mask (default: 2)
#' }
#'
#' @return A LAS object with topographic features classified into classes 19-21
#'
#' @examples
#' \dontrun{
#' library(lidR)
#' library(terra)
#' library(data.table)
#' library(spanner)
#' library(RANN)
#' 
#' # Read and normalize LAS file
#' las <- readLAS("myfile.las")
#' las <- classify_ground(las, csf())
#' dtm <- rasterize_terrain(las, res = 0.5, algorithm = tin())
#' las_norm <- normalize_height(las, dtm)
#' 
#' # Classify all topographic features
#' las_classified <- classify_topo_features(las_norm, 
#'                                           features = c("cliff", "hoodoo", "depression"))
#' 
#' # Classify only cliffs
#' las_classified <- classify_topo_features(las_norm, features = "cliff")
#' 
#' # Custom parameters using parameter lists
#' las_classified <- classify_topo_features(
#'   las_norm,
#'   features = c("cliff", "depression"),
#'   cliff_params = list(slope_thr = 60, edge_buffer_m = 10),
#'   depression_params = list(min_depth = 0.5, min_area = 5.0)
#' )
#' }
#'
#' @export
classify_topo_features <- function(
    las_norm,
    features = c("cliff", "hoodoo", "depression"),
    classes = list(cliff = 19, hoodoo = 20, depression = 21),
    dtm_res = 0.5,
    cliff_params = list(),
    hoodoo_params = list(),
    depression_params = list(),
    ncpu = 8,
    verbose = TRUE
) {
  
  # ===========================================================================
  # Default parameters for each feature type
  # ===========================================================================
  
  # Cliff detection defaults
  cliff_defaults <- list(
    slope_thr = 65,
    rescue_h = 0.35,
    mask_expand_k = 9,
    edge_buffer_m = 8,
    edge_max_h_m = 3.0,
    iso_min_h = 5.0,
    iso_search_r = 3.0,
    iso_z_tol = 2.0,
    iso_min_neighbors = 10,
    max_height = 35
  )
  
  # Hoodoo detection defaults
  hoodoo_defaults <- list(
    min_height = 6,
    top_planarity = 0.30,
    cap_break = 1.20,
    axis_wobble = 0.18,
    near_steep_m = 15
  )
  
  # Depression detection defaults
  depression_defaults <- list(
    min_depth = 0.3,
    max_depth = 5.0,
    min_area = 2.0,
    expand_cells = 2
  )
  
  # Merge user parameters with defaults
  cliff_p <- modifyList(cliff_defaults, cliff_params)
  hoodoo_p <- modifyList(hoodoo_defaults, hoodoo_params)
  dep_p <- modifyList(depression_defaults, depression_params)
  
  # Extract class assignments (with fallback)
  cliff_class <- if (!is.null(classes$cliff)) classes$cliff else 19
  hoodoo_class <- if (!is.null(classes$hoodoo)) classes$hoodoo else 20
  depression_class <- if (!is.null(classes$depression)) classes$depression else 21
  
  # ===========================================================================
  # Input validation
  # ===========================================================================
  if (!inherits(las_norm, "LAS")) {
    stop("Input must be a LAS object")
  }
  
  if (lidR::is.empty(las_norm)) {
    stop("Input LAS is empty")
  }
  
  if (!"Zref" %in% names(las_norm@data)) {
    stop("Input LAS must have a 'Zref' column (normalized LAS with ground reference)")
  }
  
  valid_features <- c("cliff", "hoodoo", "depression")
  if (!all(features %in% valid_features)) {
    invalid <- features[!features %in% valid_features]
    stop("Invalid features: ", paste(invalid, collapse = ", "), 
         ". Valid options: ", paste(valid_features, collapse = ", "))
  }
  
  if (verbose) {
    message("=== Classifying topographic features ===")
    message("Features to classify: ", paste(features, collapse = ", "))
  }
  
  # ===========================================================================
  # Helper functions
  # ===========================================================================
  
  # Extract last column from terra::extract output
  extract_last_col <- function(x) x[[ncol(x)]]
  
  # Create focal kernel from distance in meters
  kernel_from_meters <- function(r, dist_m) {
    cell_m <- terra::res(r)[1]
    k <- max(1, round(dist_m / cell_m))
    matrix(1, 2*k + 1, 2*k + 1)
  }
  
  # ===========================================================================
  # Hoodoo detection helper functions
  # ===========================================================================
  
  segment_vertical_objects <- function(las_input, min_height = 3, res = 0.5) {
    
    if (lidR::is.empty(las_input)) stop("las_input is empty")
    
    # Build CHM (pitfree with fallback)
    chm <- try(
      lidR::rasterize_canopy(
        las_input,
        res = res,
        algorithm = lidR::pitfree(thresholds = c(0, 10, 20), max_edge = c(0, 1.5))
      ),
      silent = TRUE
    )
    
    if (inherits(chm, "try-error") || is.null(chm) || all(is.na(terra::values(chm)))) {
      if (verbose) message("  pitfree CHM failed/empty; falling back to p2r() CHM")
      chm <- lidR::rasterize_canopy(las_input, res = res, algorithm = lidR::p2r())
    }
    
    if (is.null(chm) || all(is.na(terra::values(chm)))) {
      stop("CHM is empty (all NA). Check point density, filters, and CRS/extent.")
    }
    
    # Variable window function for LMF
    ws_log <- function(x, ws_min = 1.0, ws_max = 7.0) {
      x <- as.numeric(x)
      x[!is.finite(x)] <- NA_real_
      
      y <- ifelse(is.na(x), ws_min,
                  ifelse(x < 0, ws_min,
                         ifelse(x < 3, 1.35,
                                ifelse(x > 28, 5,
                                       2.51503 + 0.00901 * x^2 - 0.08297 * x^2.5))))
      
      y[!is.finite(y) | y <= 0] <- ws_min
      y <- pmin(ws_max, pmax(ws_min, y))
      y
    }
    
    # Locate tree tops from CHM
    tree_tops <- lidR::locate_trees(
      chm,
      algorithm = lidR::lmf(ws = function(z) ws_log(z, ws_min = 1.0, ws_max = 7.0),
                            hmin = min_height)
    )
    
    if (is.null(tree_tops) || nrow(tree_tops) == 0) {
      if (verbose) message("  No tree tops found; returning unsegmented LAS with treeID = NA")
      las_input@data$treeID <- NA_integer_
      return(las_input)
    }
    
    # Segment using Dalponte algorithm
    lidR::segment_trees(las_input, lidR::dalponte2016(chm, tree_tops))
  }
  
  calculate_object_metrics <- function(las_seg, eigenmetric_radius = 1.5, ncpu_local = 8) {
    
    las_eigen <- as.data.table(spanner::eigen_metrics(las_seg, radius = eigenmetric_radius, ncpu = ncpu_local))
    
    if (nrow(las_eigen) == nrow(las_seg@data)) {
      las_eigen[, treeID := las_seg@data$treeID]
      las_eigen[, `:=`(Z = las_seg@data$Z, X = las_seg@data$X, Y = las_seg@data$Y)]
    } else {
      stop("eigen_metrics() row count mismatch.")
    }
    
    # Calculate metrics for hoodoo classification
    obj_metrics <- las_eigen[!is.na(treeID), .(
      n_points = .N,
      z_range = max(Z, na.rm = TRUE) - min(Z, na.rm = TRUE),
      
      # TOP PLANARITY (hoodoo cap signature)
      top_planarity = {
        z_thresh <- quantile(Z, 0.85, na.rm = TRUE)
        top_pts <- Planarity[Z >= z_thresh]
        if (length(top_pts) > 5) mean(top_pts, na.rm = TRUE) else NA_real_
      },
      
      # AXIS WOBBLE (verticality consistency)
      axis_wobble = {
        if (.N < 30) NA_real_ else {
          z_bins <- cut(Z, breaks = 5, labels = FALSE)
          bin_vert <- tapply(Verticality, z_bins, mean, na.rm = TRUE)
          sd(bin_vert, na.rm = TRUE)
        }
      },
      
      # CAP BREAK (radial spread)
      cap_break_ratio = {
        if (.N < 30) NA_real_ else {
          x_c <- mean(X); y_c <- mean(Y)
          rad <- sqrt((X - x_c)^2 + (Y - y_c)^2)
          z_40 <- quantile(Z, 0.40, na.rm = TRUE)
          z_80 <- quantile(Z, 0.80, na.rm = TRUE)
          lower <- mean(rad[Z <= z_40], na.rm = TRUE)
          upper <- mean(rad[Z >= z_80], na.rm = TRUE)
          if (lower > 0) upper / lower else NA_real_
        }
      },
      
      x_mean = mean(X),
      y_mean = mean(Y)
    ), by = treeID]
    
    obj_metrics
  }
  
  classify_hoodoos_internal <- function(obj_metrics,
                                        min_height = 5,
                                        top_planarity_thr = 0.25,
                                        cap_break_thr = 1.15,
                                        axis_wobble_thr = 0.15) {
    
    obj_metrics[, is_hoodoo := (
      z_range >= min_height &
        !is.na(top_planarity) & top_planarity >= top_planarity_thr &
        !is.na(cap_break_ratio) & cap_break_ratio >= cap_break_thr &
        !is.na(axis_wobble) & axis_wobble >= axis_wobble_thr
    )]
    
    if (verbose) {
      message(sprintf("  Classified %d objects as hoodoos (out of %d)",
                      sum(obj_metrics$is_hoodoo), nrow(obj_metrics)))
    }
    obj_metrics
  }
  
  create_hoodoo_mask <- function(las_seg, obj_metrics, res = 1) {
    
    hoodoo_ids <- obj_metrics[is_hoodoo == TRUE, treeID]
    
    if (length(hoodoo_ids) == 0) {
      if (verbose) message("  No hoodoos detected")
      r <- terra::rast(terra::ext(las_seg), res = res, crs = terra::crs(las_seg))
      terra::values(r) <- 0
      return(r)
    }
    
    las_hoodoos <- lidR::filter_poi(las_seg, treeID %in% hoodoo_ids)
    if (lidR::is.empty(las_hoodoos)) {
      r <- terra::rast(terra::ext(las_seg), res = res, crs = terra::crs(las_seg))
      terra::values(r) <- 0
      return(r)
    }
    
    hoodoo_mask <- lidR::rasterize_canopy(las_hoodoos, res = res, algorithm = lidR::p2r())
    hoodoo_mask <- !is.na(hoodoo_mask)
    terra::focal(hoodoo_mask, w = matrix(1, 5, 5), fun = max, na.rm = TRUE)
  }
  
  # ===========================================================================
  # MAIN PROCESSING
  # ===========================================================================
  
  n_original <- nrow(las_norm@data)
  
  if (verbose) {
    message("=== Classifying topographic features ===")
    message("Features to classify: ", paste(features, collapse = ", "))
    message(sprintf("Input: %d points", n_original))
  }
  
  # Initialize classification tracking
  cliff_idx <- integer(0)
  hoodoo_idx <- integer(0)
  depression_idx <- integer(0)
  
  # ===========================================================================
  # Build shared resources (DTM, slope) if needed
  # Avoid redundant work - build once and reuse
  # ===========================================================================
  
  needs_dtm <- any(c("cliff", "depression") %in% features)
  dtm_ground <- NULL
  slope_deg <- NULL
  steep_mask <- NULL
  
  if (needs_dtm) {
    
    if (verbose) message("\n=== BUILDING SHARED RESOURCES ===")
    if (verbose) message("Building DTM from ground points...")
    
    # Get ground points (near Z=0 in normalized LAS)
    las_ground <- lidR::filter_poi(las_norm, Z <= 0.5)
    
    if (lidR::is.empty(las_ground) || nrow(las_ground@data) < 100) {
      warning("Insufficient ground points for DTM. Cannot perform terrain-based detection.")
      needs_dtm <- FALSE
    } else {
      
      # Create DTM using ground elevations (Zref)
      las_ground_elev <- las_ground
      las_ground_elev@data$Z <- las_ground_elev@data$Zref
      
      dtm_ground <- lidR::rasterize_terrain(las_ground_elev, res = dtm_res, algorithm = lidR::tin())
      
      # Build slope if needed for cliff detection
      if ("cliff" %in% features) {
        if (verbose) message("Calculating slope...")
        slope_deg <- terra::terrain(dtm_ground, v = "slope", unit = "degrees", neighbors = 8)
        
        steep_mask <- slope_deg >= cliff_p$slope_thr
        steep_mask <- terra::focal(steep_mask, w = matrix(1, cliff_p$mask_expand_k, cliff_p$mask_expand_k),
                                    fun = max, na.rm = TRUE)
      }
    }
  }
  
  # ===========================================================================
  # CLIFF DETECTION
  # ===========================================================================
  
  if ("cliff" %in% features && needs_dtm && !is.null(steep_mask)) {
    
    if (verbose) message("\n=== CLIFF DETECTION ===")
    
    # Extract steep mask values for all points
    if (verbose) message("1) Identifying steep-terrain points...")
    ex_steep <- terra::extract(steep_mask, cbind(las_norm@data$X, las_norm@data$Y))
    is_steep <- extract_last_col(ex_steep)
    is_steep[is.na(is_steep)] <- 0
    
    # Near-ground rescue on steep terrain
    idx_rescue <- which(
      is_steep == 1 &
        is.finite(las_norm@data$Z) &
        las_norm@data$Z >= -0.15 & las_norm@data$Z <= cliff_p$rescue_h
    )
    
    if (length(idx_rescue) > 0 && verbose) {
      message(sprintf("  Rescued %d near-ground points on steep terrain", length(idx_rescue)))
    }
    
    # Identify steep-mask points for cliff classification
    idx_steep <- which(is_steep == 1)
    
    # ===========================================================================
    # Cliff-edge buffer filter (STAGE 1)
    # ===========================================================================
    if (verbose) message("2) STAGE 1: Cliff-edge buffer cleanup...")
    
    w_edge <- kernel_from_meters(steep_mask, cliff_p$edge_buffer_m)
    edge_zone <- terra::focal(steep_mask, w = w_edge, fun = max, na.rm = TRUE)
    
    ex_edge <- terra::extract(edge_zone, cbind(las_norm@data$X, las_norm@data$Y))
    in_edge <- extract_last_col(ex_edge)
    in_edge[is.na(in_edge)] <- 0
    
    edge_remove_idx <- which(in_edge == 1 & las_norm@data$Z > cliff_p$edge_max_h_m)
    
    if (verbose && length(edge_remove_idx) > 0) {
      message(sprintf("  Edge-buffer filter: flagged %d points", length(edge_remove_idx)))
    }
    
    # ===========================================================================
    # Isolated high-point filter (STAGE 2)
    # ===========================================================================
    if (verbose) message("3) STAGE 2: Isolated high-point filter...")
    
    dt <- data.table::as.data.table(las_norm@data)
    dt[, ptidx := .I]
    
    candidates <- dt[Z > cliff_p$iso_min_h]
    
    iso_remove_idx <- integer(0)
    
    if (nrow(candidates) > 0) {
      
      n_cand <- nrow(candidates)
      if (verbose) message(sprintf("  Testing %d candidate high points...", n_cand))
      
      # Reference set: all points above 1 m
      ref <- dt[Z > 1.0]
      
      query_xy <- as.matrix(candidates[, .(X, Y)])
      ref_xy   <- as.matrix(ref[, .(X, Y)])
      
      K <- min(nrow(ref), max(50L, cliff_p$iso_min_neighbors * 5L))
      
      if (verbose) {
        message(sprintf("  Running kd-tree search (K=%d, %d query pts, %d ref pts)...",
                        K, n_cand, nrow(ref)))
      }
      
      nn_result <- RANN::nn2(data = ref_xy, query = query_xy, k = K,
                             searchtype = "radius", radius = cliff_p$iso_search_r)
      
      # Count height-compatible neighbors
      cand_z <- candidates$Z
      nn_idx <- nn_result$nn.idx
      nn_dist <- nn_result$nn.dists
      
      nn_idx[nn_idx == 0] <- NA_integer_
      
      ref_z_vec <- ref$Z
      nn_z <- matrix(ref_z_vec[nn_idx], nrow = n_cand, ncol = K)
      
      zdiff <- abs(nn_z - cand_z)
      valid <- (nn_dist <= cliff_p$iso_search_r) & (zdiff <= cliff_p$iso_z_tol) & (!is.na(nn_z))
      
      neighbor_counts <- rowSums(valid, na.rm = TRUE)
      
      # Subtract self
      self_in_ref <- (candidates$Z > 1.0)
      neighbor_counts <- neighbor_counts - as.integer(self_in_ref)
      neighbor_counts <- pmax(neighbor_counts, 0L)
      
      iso_remove_idx <- candidates[neighbor_counts < cliff_p$iso_min_neighbors, ptidx]
      
      if (verbose && length(iso_remove_idx) > 0) {
        message(sprintf("  Isolated-point filter: flagged %d points (out of %d tested)",
                        length(iso_remove_idx), n_cand))
      }
      
    } else {
      if (verbose) message("  No candidate high points to test")
    }
    
    # ===========================================================================
    # Final height cap for cliffs
    # ===========================================================================
    if (verbose) message("4) Final height cap...")
    
    height_cap_idx <- which(las_norm@data$Z > cliff_p$max_height)
    
    if (verbose && length(height_cap_idx) > 0) {
      message(sprintf("  Height cap: flagged %d points above %d m",
                      length(height_cap_idx), cliff_p$max_height))
    }
    
    # Combine all cliff indices
    cliff_idx <- unique(c(idx_steep, edge_remove_idx, iso_remove_idx, height_cap_idx))
    
    if (verbose) {
      message(sprintf("  Total cliff points identified: %d", length(cliff_idx)))
    }
  }
  
  # ===========================================================================
  # HOODOO DETECTION
  # ===========================================================================
  
  if ("hoodoo" %in% features) {
    
    if (!("cliff" %in% features)) {
      warning("Hoodoo detection requires 'cliff' detection to run first. Skipping hoodoos.")
    } else if (is.null(steep_mask)) {
      warning("No steep mask available. Skipping hoodoo detection.")
    } else {
      
      if (verbose) message("\n=== HOODOO DETECTION ===")
      
      if (verbose) message("1) Building near-steep zone...")
      
      w_near <- kernel_from_meters(steep_mask, hoodoo_p$near_steep_m)
      near_steep_zone <- terra::focal(steep_mask, w = w_near, fun = max, na.rm = TRUE)
      
      ex_near <- terra::extract(near_steep_zone, cbind(las_norm@data$X, las_norm@data$Y))
      in_near_steep <- extract_last_col(ex_near)
      in_near_steep[is.na(in_near_steep)] <- 0
      
      n_near <- sum(in_near_steep == 1, na.rm = TRUE)
      if (verbose) message(sprintf("  Points in near-steep zone: %d", n_near))
      
      if (n_near > 100) {
        
        las_near <- lidR::filter_poi(las_norm, in_near_steep == 1)
        
        if (!lidR::is.empty(las_near) && nrow(las_near@data) > 100) {
          
          if (verbose) message("2) Segmenting vertical objects...")
          
          las_seg <- tryCatch(
            segment_vertical_objects(las_near, min_height = hoodoo_p$min_height),
            error = function(e) {
              if (verbose) message("  Hoodoo segmentation failed: ", e$message)
              NULL
            }
          )
          
          if (!is.null(las_seg)) {
            
            if (verbose) message("3) Calculating object metrics...")
            
            obj_metrics <- tryCatch(
              calculate_object_metrics(las_seg, eigenmetric_radius = 1.5, ncpu_local = ncpu),
              error = function(e) {
                if (verbose) message("  Hoodoo metrics failed: ", e$message)
                NULL
              }
            )
            
            if (!is.null(obj_metrics) && nrow(obj_metrics) > 0) {
              
              if (verbose) message("4) Classifying hoodoos...")
              
              obj_metrics <- classify_hoodoos_internal(
                obj_metrics, hoodoo_p$min_height,
                hoodoo_p$top_planarity, hoodoo_p$cap_break, hoodoo_p$axis_wobble
              )
              
              hoodoo_ids <- obj_metrics[is_hoodoo == TRUE, treeID]
              
              if (length(hoodoo_ids) > 0) {
                
                if (verbose) message("5) Creating hoodoo mask...")
                
                hoodoo_mask <- create_hoodoo_mask(las_seg, obj_metrics, res = 1)
                hoodoo_mask_resampled <- terra::resample(hoodoo_mask, steep_mask, method = "near")
                hoodoo_final <- hoodoo_mask_resampled & near_steep_zone
                
                ex_hoodoo <- terra::extract(hoodoo_final, cbind(las_norm@data$X, las_norm@data$Y))
                is_hoodoo_pt <- extract_last_col(ex_hoodoo)
                is_hoodoo_pt[is.na(is_hoodoo_pt)] <- 0
                
                hoodoo_idx <- which(is_hoodoo_pt == 1)
                
                if (verbose) {
                  message(sprintf("  Total hoodoo points identified: %d", length(hoodoo_idx)))
                }
              } else {
                if (verbose) message("  No hoodoos classified")
              }
            }
          }
        } else {
          if (verbose) message("  Too few points in near-steep zone for hoodoo analysis")
        }
      } else {
        if (verbose) message("  Too few points in near-steep zone for hoodoo detection")
      }
    }
  }
  
  # ===========================================================================
  # DEPRESSION DETECTION
  # Based on level-set algorithm from Wu et al. (2019)
  # https://onlinelibrary.wiley.com/doi/10.1111/1752-1688.12689
  # ===========================================================================
  
  if ("depression" %in% features && needs_dtm && !is.null(dtm_ground)) {
    
    if (verbose) message("\n=== DEPRESSION DETECTION ===")
    
    # Use pre-built DTM (already created in shared resources section)
    if (verbose) message("1) Filling depressions using level-set algorithm...")
    
    # Convert to matrix for processing
    dtm_vals <- terra::values(dtm_ground, mat = TRUE)
    nr <- nrow(dtm_vals)
    nc <- ncol(dtm_vals)
    
    # Initialize filled DEM (copy of original)
    dtm_filled <- dtm_vals
    
    # Iterative depression filling using priority queue approach
    # This implements the level-set method from Wu et al. (2019)
    max_iterations <- 1000
    epsilon <- 0.001  # convergence threshold in meters
    
    for (iter in 1:max_iterations) {
      changed <- FALSE
      
      # For each cell, if it's lower than all drainable neighbors, raise it
      for (i in 2:(nr-1)) {
        for (j in 2:(nc-1)) {
          
          if (is.na(dtm_filled[i, j])) next
          
          # Get 8-neighbor values
          neighbors <- c(
            dtm_filled[i-1, j-1], dtm_filled[i-1, j], dtm_filled[i-1, j+1],
            dtm_filled[i, j-1],                       dtm_filled[i, j+1],
            dtm_filled[i+1, j-1], dtm_filled[i+1, j], dtm_filled[i+1, j+1]
          )
          
          # Remove NA neighbors
          neighbors <- neighbors[!is.na(neighbors)]
          
          if (length(neighbors) > 0) {
            # Minimum neighbor elevation
            min_neighbor <- min(neighbors)
            
            # If current cell is lower than lowest neighbor, it's in a depression
            if (dtm_filled[i, j] < min_neighbor) {
              # Raise to minimum neighbor level (level-set approach)
              dtm_filled[i, j] <- min_neighbor
              changed <- TRUE
            }
          }
        }
      }
      
      # Check convergence
      if (!changed) {
        if (verbose) message(sprintf("  Converged after %d iterations", iter))
        break
      }
      
      if (iter == max_iterations && verbose) {
        message("  Warning: Maximum iterations reached")
      }
    }
    
    # ===========================================================================
    # Calculate depression depth
    # ===========================================================================
    if (verbose) message("2) Calculating depression depth...")
    
    # Depression depth = filled - original
    dep_depth <- dtm_filled - dtm_vals
    
    # Create raster of depression depth
    dep_depth_rast <- dtm_ground
    terra::values(dep_depth_rast) <- dep_depth
    
    # ===========================================================================
    # Identify significant depressions
    # ===========================================================================
    if (verbose) message("3) Identifying significant depressions...")
    
    # Threshold by depth
    dep_mask <- (dep_depth >= dep_p$min_depth) & (dep_depth <= dep_p$max_depth)
    dep_mask[is.na(dep_mask)] <- FALSE
    
    dep_mask_rast <- dtm_ground
    terra::values(dep_mask_rast) <- as.integer(dep_mask)
    
    # Filter by minimum area using connected component labeling
    # Group connected depression cells
    if (sum(dep_mask) > 0) {
      
      # Use focal to expand depressions slightly (catch edge points)
      if (dep_p$expand_cells > 0) {
        w_dep <- matrix(1, 2*dep_p$expand_cells + 1, 2*dep_p$expand_cells + 1)
        dep_mask_rast <- terra::focal(dep_mask_rast, w = w_dep, fun = max, na.rm = TRUE)
      }
      
      # Patch/clump analysis to get connected depression areas
      dep_patches <- terra::patches(dep_mask_rast, directions = 8, zeroAsNA = TRUE)
      
      # Calculate area of each patch
      cell_area <- terra::res(dep_patches)[1]^2  # m²
      
      if (!is.null(dep_patches) && !all(is.na(terra::values(dep_patches)))) {
        
        # Get patch statistics
        patch_freq <- terra::freq(dep_patches)
        patch_freq <- patch_freq[!is.na(patch_freq$value), ]
        
        if (nrow(patch_freq) > 0) {
          patch_freq$area_m2 <- patch_freq$count * cell_area
          
          # Keep only patches above minimum area
          keep_patches <- patch_freq$value[patch_freq$area_m2 >= dep_p$min_area]
          
          if (length(keep_patches) > 0) {
            
            if (verbose) {
              message(sprintf("  Found %d depression patches (total cells: %d, mean area: %.1f m²)",
                              length(keep_patches),
                              sum(patch_freq$count[patch_freq$value %in% keep_patches]),
                              mean(patch_freq$area_m2[patch_freq$value %in% keep_patches])))
            }
            
            # Create final depression mask
            dep_final_vals <- terra::values(dep_patches, mat = FALSE)
            dep_final <- dep_final_vals %in% keep_patches
            dep_final[is.na(dep_final)] <- FALSE
            
            dep_final_rast <- dep_patches
            terra::values(dep_final_rast) <- as.integer(dep_final)
            
            # ===========================================================================
            # Extract points in depressions
            # ===========================================================================
            if (verbose) message("4) Extracting points in depression areas...")
            
            ex_dep <- terra::extract(dep_final_rast, cbind(las_norm@data$X, las_norm@data$Y))
            in_dep <- extract_last_col(ex_dep)
            in_dep[is.na(in_dep)] <- 0
            
            depression_idx <- which(in_dep == 1)
            
            if (verbose) {
              message(sprintf("  Total depression points identified: %d", length(depression_idx)))
            }
            
          } else {
            if (verbose) message("  No depressions meet minimum area threshold")
            depression_idx <- integer(0)
          }
        } else {
          if (verbose) message("  No depression patches found")
          depression_idx <- integer(0)
        }
      } else {
        if (verbose) message("  No depression patches detected")
        depression_idx <- integer(0)
      }
      
    } else {
      if (verbose) message("  No depressions meet depth threshold")
      depression_idx <- integer(0)
    }
  }
  
  # ===========================================================================
  # APPLY CLASSIFICATIONS
  # ===========================================================================
  
  if (verbose) message("\n=== APPLYING CLASSIFICATIONS ===")
  
  # Remove overlaps: hoodoo takes precedence over cliff
  # (hoodoos are a more specific feature)
  cliff_only_idx <- setdiff(cliff_idx, hoodoo_idx)
  
  if (length(cliff_only_idx) > 0) {
    las_norm@data$Classification[cliff_only_idx] <- cliff_class
    if (verbose) {
      message(sprintf("  Classified %d points as cliff (class %d)", 
                      length(cliff_only_idx), cliff_class))
    }
  }
  
  if (length(hoodoo_idx) > 0) {
    las_norm@data$Classification[hoodoo_idx] <- hoodoo_class
    if (verbose) {
      message(sprintf("  Classified %d points as hoodoo (class %d)", 
                      length(hoodoo_idx), hoodoo_class))
    }
  }
  
  if (length(depression_idx) > 0) {
    las_norm@data$Classification[depression_idx] <- depression_class
    if (verbose) {
      message(sprintf("  Classified %d points as depression (class %d)", 
                      length(depression_idx), depression_class))
    }
  }
  
  # ===========================================================================
  # SUMMARY
  # ===========================================================================
  
  n_classified <- length(cliff_only_idx) + length(hoodoo_idx) + length(depression_idx)
  
  if (verbose) {
    message("\n=== SUMMARY ===")
    message(sprintf("Total points processed:      %d", n_original))
    message(sprintf("Cliff points (class %d):      %d", cliff_class, length(cliff_only_idx)))
    message(sprintf("Hoodoo points (class %d):     %d", hoodoo_class, length(hoodoo_idx)))
    message(sprintf("Depression points (class %d): %d", depression_class, length(depression_idx)))
    message(sprintf("Total classified:            %d (%.1f%%)", 
                    n_classified, 100 * n_classified / n_original))
  }
  
  return(las_norm)
}


###############################################################################
# classify_topo_as_ground
#
# Reclassifies topographic feature classes (19, 20, 21) as ground (class 2)
###############################################################################

#' Reclassify Topographic Features as Ground
#'
#' Converts cliff (19), hoodoo (20), and depression (21) points to ground 
#' class (2). This is useful for treating topographic features as terrain 
#' rather than vegetation.
#'
#' @param las LAS object with topographic features classified
#' @param classes Numeric vector of classes to reclassify (default: c(19, 20, 21))
#' @param ground_class Numeric class to assign (default: 2 for ground)
#' @param verbose Print messages (default: TRUE)
#'
#' @return LAS object with reclassified points
#'
#' @examples
#' \dontrun{
#' # Reclassify all topographic features as ground
#' las <- classify_topo_as_ground(las)
#' 
#' # Reclassify only cliffs and hoodoos
#' las <- classify_topo_as_ground(las, classes = c(19, 20))
#' }
#'
#' @export
classify_topo_as_ground <- function(las, 
                                     classes = c(19, 20, 21), 
                                     ground_class = 2,
                                     verbose = TRUE) {
  
  # Count points in each class before reclassification
  idx_to_reclassify <- which(las@data$Classification %in% classes)
  n_reclassify <- length(idx_to_reclassify)
  
  if (n_reclassify == 0) {
    if (verbose) {
      message("No topographic feature points found to reclassify")
    }
    return(las)
  }
  
  # Get breakdown by class
  if (verbose) {
    message("\n=== Reclassifying Topographic Features as Ground ===")
    for (cls in classes) {
      n_cls <- sum(las@data$Classification == cls)
      if (n_cls > 0) {
        message(sprintf("  Class %d: %d points -> ground (class %d)", 
                        cls, n_cls, ground_class))
      }
    }
  }
  
  # Reclassify
  las@data$Classification[idx_to_reclassify] <- ground_class
  
  if (verbose) {
    message(sprintf("Total reclassified: %d points", n_reclassify))
  }
  
  return(las)
}


###############################################################################
# classify_voxel_neighborhood
#
# Voxel-based classification smoothing using majority neighbor voting
###############################################################################

#' Voxel-based Classification Smoothing
#'
#' Reclassifies isolated points based on the majority class of neighboring 
#' voxels. Similar to lidR::ivf() but for classification smoothing. For each 
#' point above a height threshold, examines neighboring voxels and reclassifies 
#' if the point has sufficient neighbors of a different class.
#'
#' @param las LAS object with Classification field
#' @param res Numeric. Voxel resolution in meters (default: 0.5)
#' @param min_neighbors Integer. Minimum number of neighbors required to 
#'   trigger reclassification (default: 5)
#' @param min_height Numeric. Minimum height (Z) to apply filter. Points 
#'   below this height are not reclassified (default: 0.5)
#' @param exclude_classes Numeric vector of classes to exclude from 
#'   reclassification (default: c(2, 9) for ground and water)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @details
#' This function implements a voxel-based classification smoothing algorithm:
#' 1. Divides the point cloud into voxels of size \code{res}
#' 2. For each point above \code{min_height}
#' 3. Finds all neighboring voxels (26-neighborhood in 3D)
#' 4. Counts points by class in neighboring voxels
#' 5. If the point's class differs from the majority neighbor class AND
#'    there are at least \code{min_neighbors} neighbors, reclassifies the 
#'    point to the majority class
#'
#' This helps remove isolated misclassified points and smooth classification 
#' boundaries.
#'
#' @return LAS object with smoothed classifications
#'
#' @examples
#' \dontrun{
#' # Basic smoothing with default parameters
#' las <- classify_voxel_neighborhood(las)
#' 
#' # More aggressive smoothing with larger voxels
#' las <- classify_voxel_neighborhood(las, res = 1.0, min_neighbors = 3)
#' 
#' # Only smooth vegetation points, exclude ground
#' las <- classify_voxel_neighborhood(las, exclude_classes = c(2))
#' }
#'
#' @export
classify_voxel_neighborhood <- function(las,
                                         res = 0.5,
                                         min_neighbors = 5,
                                         min_height = 0.5,
                                         exclude_classes = c(2, 9),
                                         verbose = TRUE) {
  
  if (verbose) {
    message("\n=== Voxel-based Classification Smoothing ===")
    message(sprintf("Voxel resolution: %.2f m", res))
    message(sprintf("Minimum neighbors: %d", min_neighbors))
    message(sprintf("Minimum height: %.2f m", min_height))
  }
  
  # Extract coordinates and classification
  xyz <- las@data[, c("X", "Y", "Z")]
  classification <- las@data$Classification
  n_original <- nrow(xyz)
  
  # Filter points above height threshold and not in excluded classes
  apply_filter <- (xyz$Z > min_height) & !(classification %in% exclude_classes)
  candidate_idx <- which(apply_filter)
  
  if (length(candidate_idx) == 0) {
    if (verbose) message("No points meet criteria for smoothing")
    return(las)
  }
  
  if (verbose) {
    message(sprintf("Processing %d points (%.1f%% of total)", 
                    length(candidate_idx), 
                    100 * length(candidate_idx) / n_original))
  }
  
  # Create voxel indices for all points
  voxel_x <- floor(xyz$X / res)
  voxel_y <- floor(xyz$Y / res)
  voxel_z <- floor(xyz$Z / res)
  
  # Create voxel key (unique identifier for each voxel)
  # Using a formula to create unique integer keys
  x_range <- max(voxel_x) - min(voxel_x) + 1
  y_range <- max(voxel_y) - min(voxel_y) + 1
  
  voxel_key <- (voxel_x - min(voxel_x)) + 
               (voxel_y - min(voxel_y)) * x_range + 
               (voxel_z - min(voxel_z)) * x_range * y_range
  
  # Build voxel lookup table
  voxel_data <- data.table::data.table(
    idx = 1:n_original,
    voxel_x = voxel_x,
    voxel_y = voxel_y,
    voxel_z = voxel_z,
    voxel_key = voxel_key,
    class = classification
  )
  
  # Create index by voxel for fast lookup
  data.table::setkey(voxel_data, voxel_key)
  
  # Track reclassifications
  n_reclassified <- 0
  reclassifications <- integer(0)
  new_classes <- integer(0)
  
  # Process each candidate point
  for (i in candidate_idx) {
    
    # Get voxel coordinates of current point
    vx <- voxel_x[i]
    vy <- voxel_y[i]
    vz <- voxel_z[i]
    current_class <- classification[i]
    
    # 26-neighborhood in 3D (all adjacent voxels)
    neighbor_classes <- integer(0)
    
    for (dx in -1:1) {
      for (dy in -1:1) {
        for (dz in -1:1) {
          # Skip center voxel
          if (dx == 0 && dy == 0 && dz == 0) next
          
          # Calculate neighbor voxel coordinates
          nvx <- vx + dx
          nvy <- vy + dy
          nvz <- vz + dz
          
          # Calculate neighbor voxel key
          nkey <- (nvx - min(voxel_x)) + 
                  (nvy - min(voxel_y)) * x_range + 
                  (nvz - min(voxel_z)) * x_range * y_range
          
          # Get points in this neighbor voxel
          neighbor_pts <- voxel_data[voxel_key == nkey]
          
          if (nrow(neighbor_pts) > 0) {
            neighbor_classes <- c(neighbor_classes, neighbor_pts$class)
          }
        }
      }
    }
    
    # If we have enough neighbors
    if (length(neighbor_classes) >= min_neighbors) {
      
      # Find majority class among neighbors
      class_counts <- table(neighbor_classes)
      majority_class <- as.integer(names(class_counts)[which.max(class_counts)])
      
      # If current class differs from majority, reclassify
      if (current_class != majority_class) {
        classification[i] <- majority_class
        n_reclassified <- n_reclassified + 1
        reclassifications <- c(reclassifications, i)
        new_classes <- c(new_classes, majority_class)
      }
    }
  }
  
  # Apply reclassifications
  if (n_reclassified > 0) {
    las@data$Classification <- classification
  }
  
  if (verbose) {
    message(sprintf("Reclassified %d points (%.1f%% of candidates)", 
                    n_reclassified, 
                    100 * n_reclassified / length(candidate_idx)))
    
    if (n_reclassified > 0) {
      # Show breakdown of reclassifications
      class_changes <- table(new_classes)
      message("Reclassification breakdown:")
      for (cls in sort(unique(new_classes))) {
        n_cls <- sum(new_classes == cls)
        message(sprintf("  -> Class %d: %d points", cls, n_cls))
      }
    }
  }
  
  return(las)
}
