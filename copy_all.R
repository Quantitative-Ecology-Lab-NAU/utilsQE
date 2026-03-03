#' Copy all files and folders from one directory to another
#'
#' This function copies all files and subdirectories from a source 
#' directory to a destination directory. The directory structure 
#' is preserved, and the destination directory is created if it 
#' does not already exist.
#'
#' @param from Character. Path to the source directory
#' @param to Character. Path to the destination directory
#' @param overwrite Logical. Overwrite existing files in destination 
#'   (default: FALSE)
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return Logical vector indicating which files were successfully copied
#'
#' @details
#' This function performs the following operations:
#' 1. Validates that the source directory exists
#' 2. Creates the destination directory if necessary
#' 3. Lists all top-level files and folders in the source
#' 4. Recursively copies all content while preserving structure
#'
#' Note: For very large directory trees, system-level tools such as 
#' `rsync` or `cp -r` may provide better performance.
#'
#' @examples
#' \dontrun{
#' # Copy all contents
#' copy_directory(
#'   from = "/home/rstudio/data",
#'   to   = "/home/rstudio/backup"
#' )
#'
#' # Copy and overwrite existing files
#' copy_directory(
#'   from = "input_folder",
#'   to   = "output_folder",
#'   overwrite = TRUE
#' )
#' }
#'
#' @export
copy_directory <- function(from,
                           to,
                           overwrite = FALSE,
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
  
  # Create destination directory if needed
  if (!dir.exists(to)) {
    if (verbose) message("Creating destination directory: ", to)
    dir.create(to, recursive = TRUE)
  }
  
  # List top-level items
  items <- list.files(from, full.names = TRUE, recursive = FALSE)
  
  if (length(items) == 0) {
    if (verbose) message("Source directory is empty.")
    return(invisible(logical(0)))
  }
  
  if (verbose) {
    message("Copying ", length(items), " item(s) from:")
    message("  ", normalizePath(from))
    message("To:")
    message("  ", normalizePath(to))
  }
  
  # Perform recursive copy
  result <- file.copy(
    from = items,
    to = to,
    recursive = TRUE,
    overwrite = overwrite
  )
  
  if (verbose) {
    message("Copy complete.")
  }
  
  return(invisible(result))
}
