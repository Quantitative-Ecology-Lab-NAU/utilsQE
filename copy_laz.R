#' Copy LAS/LAZ files from one directory tree to another with a progress bar
#'
#' This function recursively finds all .las and .laz files under a source
#' directory and copies them into a destination directory, preserving the
#' relative folder structure. A progress bar is shown during copying.
#'
#' @param from Character. Path to the source directory
#' @param to Character. Path to the destination directory
#' @param overwrite Logical. Overwrite existing files in destination
#'   (default: FALSE)
#' @param recursive Logical. Search subdirectories (default: TRUE)
#' @param verbose Logical. Print summary messages (default: TRUE)
#'
#' @return Invisibly returns a data.frame with source, destination, and success flag
#'
#' @details
#' This function performs the following operations:
#' 1. Validates that the source directory exists
#' 2. Creates the destination directory if necessary
#' 3. Recursively lists LAS/LAZ files under `from`
#' 4. Computes each file's destination path by preserving relative structure
#' 5. Copies files one-by-one, creating destination subfolders as needed
#' 6. Displays a progress bar and returns per-file copy status
#'
#' @examples
#' \dontrun{
#' copy_laz(
#'   from = "/data-store/home/bi0m3trics/lidar/source_tiles",
#'   to   = "/data-store/home/bi0m3trics/lidar/laz_subset",
#'   overwrite = FALSE
#' )
#' }
#'
#' @export
copy_laz <- function(from,
                     to,
                     overwrite = FALSE,
                     recursive = TRUE,
                     verbose = TRUE) {
  
  # Validate inputs
  if (!is.character(from) || length(from) != 1) {
    stop("'from' must be a single character string path", call. = FALSE)
  }
  
  if (!is.character(to) || length(to) != 1) {
    stop("'to' must be a single character string path", call. = FALSE)
  }
  
  if (!dir.exists(from)) {
    stop("Source directory does not exist: ", from, call. = FALSE)
  }
  
  if (!dir.exists(to)) {
    if (verbose) message("Creating destination directory: ", to)
    dir.create(to, recursive = TRUE, showWarnings = FALSE)
  }
  
  from_norm <- normalizePath(from, winslash = "/", mustWork = TRUE)
  to_norm   <- normalizePath(to, winslash = "/", mustWork = TRUE)
  
  # List LAS/LAZ files
  files <- list.files(
    from_norm,
    pattern = "\\.(las|laz)$",
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )
  
  if (length(files) == 0) {
    if (verbose) message("No .las/.laz files found under: ", from_norm)
    out <- data.frame(from = character(0), to = character(0), ok = logical(0))
    return(invisible(out))
  }
  
  # Compute relative paths and destination paths
  rel <- sub(paste0("^", from_norm, "/?"), "", normalizePath(files, winslash = "/"))
  dest <- file.path(to_norm, rel)
  
  if (verbose) {
    message("Found ", length(files), " LAS/LAZ file(s). Copying to:")
    message("  ", to_norm)
  }
  
  pb <- utils::txtProgressBar(min = 0, max = length(files), style = 3)
  on.exit(close(pb), add = TRUE)
  
  ok <- logical(length(files))
  
  for (i in seq_along(files)) {
    # Ensure destination subfolder exists
    dir.create(dirname(dest[i]), recursive = TRUE, showWarnings = FALSE)
    
    ok[i] <- file.copy(
      from = files[i],
      to = dest[i],
      overwrite = overwrite
    )
    
    utils::setTxtProgressBar(pb, i)
  }
  
  if (verbose) {
    message("\nCopy complete: ", sum(ok), "/", length(ok), " succeeded.")
    if (any(!ok)) {
      message("Failures: ", sum(!ok))
    }
  }
  
  out <- data.frame(from = files, to = dest, ok = ok, stringsAsFactors = FALSE)
  return(invisible(out))
}
