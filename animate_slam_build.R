# ============================================================================
# Animated SLAM Point Cloud Build — GeoSLAM ZEB Horizon
# lidR + rgl + magick
#
# Shows the point cloud being "built" one second at a time:
#   • Current second's points → white
#   • All previously accumulated points → colored by a user-chosen variable
#   • Scanner trajectory path → red points
# Saved as an animated GIF.
# ============================================================================

# ---- USER CONFIGURATION ---------------------------------------------------

laz_file   <- "c:/Users/ajsm/Desktop/SLAM_vis/2021-10-15_13-22-53.laz"
traj_file  <- "c:/Users/ajsm/Desktop/SLAM_vis/2021-10-15_13-22-53.txt"
output_gif <- "c:/Users/ajsm/Desktop/SLAM_vis/slam_build.gif"

# Color mapping
color_var  <- "Z"                    # attribute to map colour to
color_pal  <- lidR::height.colors    # palette function  (takes n → n colours)
n_colors   <- 50                     # number of bins in the palette

# Scanner path
path_color <- "red"                  # colour for scanner trajectory points
path_size  <- 4                      # point size for scanner path

# Animation / display
gif_fps    <- 10                     # playback speed (frames per second)
pt_size    <- 1                      # rgl point size
bg_color   <- "black"                # background colour
new_color  <- "white"                # colour for the "active" second's points
show_label <- TRUE                   # overlay elapsed-time label on each frame

# RGL display size (pixels) — controls the rendered scene dimensions
display_w  <- 800
display_h  <- 800

# View — custom set_view() function.  Named viewpoints:
#   "top"  (top-down)       "bottom" (bottom-up)
#   "front"                 "back"
#   "left"                  "right"
#   "front isometric"       "back isometric"
#   "custom"  → supply a 4×4 userMatrix via view_custom_um
view_point     <- "top"
view_FOV       <- 0                        # 0 = orthographic; >0 = perspective
view_zoom      <- 0.75
view_custom_um <- NULL                     # only used when view_point = "custom"
                                           # e.g. rotationMatrix(pi/6, 1, 0.5, 0)

# Performance — set NULL to keep ALL points (will be very slow for large clouds)
max_pts_per_sec <- 5000

# ---- LIBRARIES -------------------------------------------------------------

library(lidR)
library(rgl)
library(data.table)
library(magick)

# ---- CUSTOM VIEW HELPER ------------------------------------------------------
#
# set_view()  — applies a named camera orientation + FOV + zoom to the current
#               rgl device, replacing the nat::nview3d() dependency.
#
#   Coordinate convention (GeoSLAM ZEB with XYZ as easting/northing/up):
#     X → right,  Y → forward (into screen),  Z → up
#
#   rgl's default camera looks along –Z in *eye* space.  The 4×4 userMatrix
#   rotates the *scene* before the camera sees it, so we build rotation
#   matrices that reorient the scene to place the desired axis toward the camera.
#
set_view <- function(viewpoint = "top", FOV = 0, zoom = 0.75,
                     custom_um = NULL) {

  # Rotation helper — 4×4 rotation about an arbitrary axis (angle in radians)
  R4 <- function(angle, x, y, z) {
    rgl::rotationMatrix(angle, x, y, z)
  }

  um <- switch(tolower(viewpoint),
    # Top-down: camera looks down –Z.  Scene identity (Z already points up).
    "top" = diag(4),

    # Bottom-up: flip scene 180° around X so Z points down toward camera.
    "bottom" = R4(pi, 1, 0, 0),

    # Front: camera faces the –Y side.  Tilt scene –90° around X.
    "front" = R4(-pi/2, 1, 0, 0),

    # Back: camera faces the +Y side.  Tilt scene +90° around X.
    "back" = R4(pi/2, 1, 0, 0),

    # Left: view from –X.  Tilt –90° around X, then +90° around Y.
    "left" = R4(pi/2, 0, 1, 0) %*% R4(-pi/2, 1, 0, 0),

    # Right: view from +X.  Tilt –90° around X, then –90° around Y.
    "right" = R4(-pi/2, 0, 1, 0) %*% R4(-pi/2, 1, 0, 0),

    # Front isometric: slight tilt to show depth — front-and-above.
    "front isometric" = R4(-pi/6, 1, 0, 0) %*% R4(pi/6, 0, 0, 1),

    # Back isometric: same idea from the opposite side.
    "back isometric" = R4(pi/6, 1, 0, 0) %*% R4(pi + pi/6, 0, 0, 1),

    # Custom: user supplies a 4×4 matrix directly.
    "custom" = {
      if (is.null(custom_um)) stop("view_point='custom' requires view_custom_um")
      custom_um
    },

    stop("Unknown viewpoint: '", viewpoint,
         "'.  Choose top/bottom/front/back/left/right/",
         "front isometric/back isometric/custom.")
  )

  par3d(FOV = FOV, zoom = zoom, userMatrix = um)
}

# ---- READ DATA --------------------------------------------------------------

cat("Reading LAZ file …\n")
las <- readLAS(laz_file)
dt  <- as.data.table(las@data)
rm(las); invisible(gc())

# Bin every point into its integer-second bucket
dt[, sec := floor(gpstime)]
unique_secs <- sort(unique(dt$sec))
n_frames    <- length(unique_secs)

cat(sprintf("  %s points  |  %d one-second bins\n",
            format(nrow(dt), big.mark = ","), n_frames))

# ---- READ TRAJECTORY ---------------------------------------------------------

cat("Reading trajectory file …\n")
# Columns: world_time x y z q0 q1 q2 q3 r g b nx ny nz roll pitch yaw
traj <- fread(traj_file, skip = 1, header = FALSE)
setnames(traj, c("world_time", "x", "y", "z",
                 "q0", "q1", "q2", "q3",
                 "r", "g", "b",
                 "nx", "ny", "nz",
                 "roll", "pitch", "yaw"))
traj[, sec := floor(world_time)]
cat(sprintf("  %s trajectory poses\n", format(nrow(traj), big.mark = ",")))

# ---- OPTIONAL SUBSAMPLE -----------------------------------------------------

if (!is.null(max_pts_per_sec)) {
  cat(sprintf("Subsampling to <= %s pts per second … ",
              format(max_pts_per_sec, big.mark = ",")))
  set.seed(42)
  dt <- dt[, {
    idx <- if (.N > max_pts_per_sec) sample(.N, max_pts_per_sec) else seq_len(.N)
    .SD[idx]
  }, by = sec]
  setkey(dt, sec)
  cat(sprintf("%s points retained\n", format(nrow(dt), big.mark = ",")))
}

# ---- PRE-COMPUTE PALETTE COLOURS --------------------------------------------

pal        <- color_pal(n_colors)
vals       <- dt[[color_var]]
vr         <- range(vals, na.rm = TRUE)
scaled     <- (vals - vr[1]) / diff(vr)
scaled     <- pmin(pmax(scaled, 0), 1)
dt[, pt_color := pal[floor(scaled * (n_colors - 1)) + 1]]

# ---- SPATIAL EXTENT (fixed across all frames) --------------------------------

xlim <- range(dt$X);  ylim <- range(dt$Y);  zlim <- range(dt$Z)

# ---- TEMP DIRECTORY FOR FRAME PNGs ------------------------------------------

tmp_dir <- file.path(tempdir(), "slam_frames")
if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
dir.create(tmp_dir, recursive = TRUE)

# ---- OPEN RGL DEVICE --------------------------------------------------------
#
#   NOTE: the rgl window must stay VISIBLE and UNOBSCURED while frames render
#   (snapshot3d reads pixels from the OpenGL buffer).
#
#   Top-down orthographic view:
#     userMatrix = diag(4)  →  identity rotation  →  camera looks along –Z
#     FOV = 0               →  orthographic (no perspective distortion)
#
#   To switch to a 3-D perspective view, change FOV and/or userMatrix, e.g.:
#     par3d(FOV = 30, userMatrix = rotationMatrix(pi/6, 1, 0.5, 0))
#

open3d(windowRect = c(50, 50, 50 + display_w, 50 + display_h), silent = TRUE)
bg3d(bg_color)

# Invisible anchor points at the 8 bounding-box corners lock the view extent
# so the camera does not re-frame as new points appear.
bbox <- as.matrix(expand.grid(X = xlim, Y = ylim, Z = zlim))
points3d(bbox[, 1], bbox[, 2], bbox[, 3],
         color = bg_color, size = 0.001)

# Apply the chosen camera view
set_view(viewpoint  = view_point,
         FOV        = view_FOV,
         zoom       = view_zoom,
         custom_um  = view_custom_um)
Sys.sleep(0.5)                       # let the window initialise on screen

# ---- RENDER LOOP -------------------------------------------------------------

cat(sprintf("Rendering %d frames …\n", n_frames))

live_pts_id <- NULL
t0          <- unique_secs[1]

for (i in seq_len(n_frames)) {

  par3d(skipRedraw = TRUE)            # batch scene changes, redraw once

  cur_sec <- unique_secs[i]
  cur     <- dt[sec == cur_sec]

  # Swap: previous white points → their palette colours
  if (!is.null(live_pts_id)) {
    try(pop3d(id = live_pts_id), silent = TRUE)
    prev <- dt[sec == unique_secs[i - 1L]]
    points3d(prev$X, prev$Y, prev$Z,
             color = prev$pt_color, size = pt_size)
  }

  # Current second → white
  live_pts_id <- points3d(cur$X, cur$Y, cur$Z,
                          color = new_color, size = pt_size)

  # Scanner trajectory: add this second's new poses (accumulates across frames)
  path_cur <- traj[sec == cur_sec]
  if (nrow(path_cur) > 0) {
    points3d(path_cur[["x"]], path_cur[["y"]], path_cur[["z"]],
             color = path_color, size = path_size)
  }

  par3d(skipRedraw = FALSE)
  Sys.sleep(0.05)                     # allow the redraw to land

  frame_path <- file.path(tmp_dir, sprintf("frame_%04d.png", i))
  snapshot3d(frame_path, fmt = "png", webshot = FALSE)

  cat(sprintf("\r  [%3d / %d]  %.0f%%", i, n_frames, 100 * i / n_frames))
}
cat("\n")
close3d()

# ---- ASSEMBLE GIF ------------------------------------------------------------

cat("Assembling GIF …\n")

frame_files <- sort(list.files(tmp_dir, "\\.png$", full.names = TRUE))
imgs <- image_read(frame_files)

# Optional elapsed-time overlay
if (show_label) {
  t_end <- max(unique_secs) - t0
  annotated <- lapply(seq_along(imgs), function(j) {
    elapsed <- unique_secs[j] - t0
    label   <- sprintf("T+%03ds / %ds", elapsed, t_end)
    image_annotate(imgs[j], label,
                   size = 22, color = "white",
                   location = "+10+10", font = "mono")
  })
  imgs <- do.call(c, annotated)
}

anim <- image_animate(imgs, fps = gif_fps, dispose = "previous")
image_write(anim, path = output_gif)

# ---- CLEANUP -----------------------------------------------------------------

unlink(tmp_dir, recursive = TRUE)
cat(sprintf("\nDone!  GIF saved to:\n  %s\n", output_gif))
