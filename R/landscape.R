# Package dependencies used below: R6, dplyr, sf, igraph, imager.
# (Only `sf`/`igraph` objects trigger the corresponding `inherits()` checks;
# the packages themselves don't need to be attached for class definition.)

#' @title landscape R6 Class
#'
#' @description
#' An R6 class representing a landscape as a set of sites, optionally paired
#' with feature data about those sites. Locations (sites) in `landscape`
#' objects can be stored as raw coordinates, spatial geometries, or
#' a neighbour graph. Optional features can be included in the form
#' of tabular data, rasters, or impassable boundaries.
#'
#' @details
#' A `landscape` stores:
#' \itemize{
#'   \item{\code{n}: the total number of stored sites.}
#'   \item{\code{sites}: site locations, as ONE of:
#'     \itemize{
#'       \item{a 2-column numeric matrix of x/y coordinates,}
#'       \item{an \code{sfc} object of \code{POINT} or \code{POLYGON} geometries, or}
#'       \item{an adjacency graph (\code{igraph} object) describing which
#'             sites neighbour which.}
#'     }}
#'   \item{\code{features}: (optional) site-level data, as ONE of:
#'     \itemize{
#'       \item{a tibble/data frame with \code{n} rows, one column per property,}
#'       \item{a list of raster objects (one per feature layer), or}
#'       \item{an \code{sf} object of boundary geometries impassable to movement.}
#'     }}
#' }
#'
#' @export
landscape <- R6::R6Class(
  classname = "landscape",

  public = list(

    #' @field n Integer. The total number of stored sites.
    n = NULL,

    #' @field sites Site locations: a 2-column x/y `matrix`, an `sfc` object
    #'   (POINT or POLYGON geometries), or an adjacency graph (`igraph`).
    sites = NULL,

    #' @field sites_type Character scalar describing how `sites` is stored.
    #'   One of `"matrix"`, `"sfc"`, or `"graph"`. Derived automatically;
    #'   not set directly.
    sites_type = NULL,

    #' @field features (Optional) site-level data: a tibble/data frame, a
    #'   list of rasters, or an `sf` object of impassable boundaries.
    features = NULL,

    #' @field features_type Character scalar describing how `features` is
    #'   stored, or `NULL` if no features are present. One of `"tibble"`,
    #'   `"raster_list"`, or `"boundaries"`. Derived automatically; not set
    #'   directly.
    features_type = NULL,

    #' @description
    #' Create a new `landscape` object.
    #'
    #' @param sites ONE of: a 2-column numeric `matrix` of x/y coordinates, an
    #'   `sfc` object of POINT or POLYGON geometries, or an adjacency graph
    #'   (`igraph` object) describing site neighbours.
    #' @param features Optional. ONE of: a tibble/data frame with one row
    #'   per site, a list of raster objects, or an `sf` object of boundary
    #'   geometries impassable to movement.
    #'
    #' @return A new `landscape` object, invisibly.
    initialize = function(sites, features = NULL) {
      private$set_sites(sites)
      private$set_features(features)
      invisible(self)
    },

    #' @description
    #' Print a summary of the landscape.
    #' @param ... Unused; included for S3 `print()` compatibility.
    print = function(...) {
      cat("<landscape>\n")
      cat("  n sites:      ", self$n, "\n", sep = "")
      cat("  sites stored as:   ", self$sites_type, "\n", sep = "")
      if (!is.null(self$features_type)) {
        cat("  features stored as:", self$features_type, "\n", sep = " ")
      } else {
        cat("  features:           none\n")
      }
      invisible(self)
    },

    #' @description
    #' Summarise the landscape's sites and features.
    #'
    #' Rough mock-up: prints a short summary and returns it invisibly as a
    #' list. The per-`features_type` branch is sketched but not fleshed
    #' out (e.g. `summary()` on a tibble is a placeholder for whatever
    #' subset of columns / stats actually matter downstream).
    #'
    #' @param ... Unused; included for S3 `summary()` compatibility.
    #' @return Invisibly, a list of summary values (also printed).
    summary = function(...) {
      cat("landscape summary\n")
      cat("------------------\n")
      cat("n sites:   ", self$n, "\n", sep = "")
      cat("sites type:", self$sites_type, "\n\n")

      feature_summary <- NULL
      if (is.null(self$features)) {
        cat("features: none\n")
      } else {
        cat("features type:", self$features_type, "\n")
        # TODO: decide what's actually useful per representation; this is
        # just a sketch of the dispatch.
        feature_summary <- switch(
          self$features_type,
          tibble      = summary(self$features),
          raster_list = paste0(length(self$features), " raster layer(s)"),
          boundaries  = paste0(nrow(self$features), " boundary feature(s)")
        )
        print(feature_summary)
      }

      invisible(list(
        n = self$n,
        sites_type = self$sites_type,
        features_type = self$features_type,
        feature_summary = feature_summary
      ))
    },

    #' @description
    #' Plot the landscape's sites.
    #'
    #' Rough mock-up: dispatches on `sites_type` to a sensible base/sf/
    #' igraph plot call. Not yet doing anything with `features` (e.g.
    #' colouring points by an intensity column when `features_type ==
    #' "tibble"`) — that's the obvious next step once the shape of
    #' `features` is nailed down.
    #'
    #' @param ... Passed on to the underlying plotting function.
    #' @return Invisibly, `self`.
    plot = function(...) {
      switch(
        self$sites_type,
        matrix = {
          plot(self$sites, pch = 16, asp = 1, xlab = "x", ylab = "y", ...)
        },
        sfc = {
          plot(self$sites, ...)
        },
        graph = {
          igraph::plot.igraph(self$sites, ...)
        }
      )
      invisible(self)
    },

    #' @description
    #' Downscale the landscape (e.g. reduce site resolution or count).
    #'
    #' @param ... TODO: define parameters (e.g. a `factor`/`scale` arg,
    #'   an aggregation function for `features`).
    #' @return TODO.
    downscale = function(...) {
      # TODO: implement
    },

    #' @description
    #' Compute inter-site distances across the landscape.
    #'
    #' @param ... TODO: define parameters (e.g. a distance metric, whether
    #'   `features`/boundaries factor into cost distance).
    #' @return TODO.
    em_dist = function(...) {
      # TODO: implement
    }
  ),

  private = list(

    # Validate `sites`, and derive `n` / `sites_type` from it.
    set_sites = function(sites) {
      if (is.matrix(sites) && is.numeric(sites)) {
        if (ncol(sites) != 2) {
          stop("`sites` matrix must have exactly 2 columns (x and y).", call. = FALSE)
        }
        self$sites_type <- "matrix"
        self$n <- nrow(sites)

      } else if (inherits(sites, "sfc")) {
        geom_types <- unique(as.character(sf::st_geometry_type(sites)))
        if (!all(geom_types %in% c("POINT", "POLYGON"))) {
          stop(
            "`sites` sfc object must contain only POINT or POLYGON geometries, ",
            "found: ", paste(geom_types, collapse = ", "), ".",
            call. = FALSE
          )
        }
        self$sites_type <- "sfc"
        self$n <- length(sites)

      } else if (inherits(sites, "igraph")) {
        self$sites_type <- "graph"
        self$n <- igraph::vcount(sites)

      } else {
        stop(
          "`sites` must be one of: a 2-column numeric matrix, an sfc object ",
          "(POINT/POLYGON geometries), or an adjacency graph (igraph object). ",
          "Instead, got an object of class: ",
          paste(class(sites), collapse = "/"), ".",
          call. = FALSE
        )
      }

      self$sites <- sites
    },

    # Validate `features` (if supplied) against `self$n`, and derive
    # `features_type`.
    set_features = function(features) {
      if (is.null(features)) {
        self$features <- NULL
        self$features_type <- NULL
        return(invisible(NULL))
      }

      if (inherits(features, "sf")) {
        self$features_type <- "boundaries"

      } else if (inherits(features, "data.frame")) { # tibbles inherit data.frame
        if (nrow(features) != self$n) {
          stop(
            "`features` must have exactly `n` (", self$n, ") rows, one per site; ",
            "got ", nrow(features), ".",
            call. = FALSE
          )
        }
        self$features_type <- "tibble"

      } else if (is.list(features) &&
                 length(features) > 0 &&
                 all(vapply(features, private$is_raster, logical(1)))) {
        self$features_type <- "raster_list"

      } else {
        stop(
          "`features` must be one of: a tibble/data frame with `n` rows, a ",
          "non-empty list of raster objects, or an sf object of boundary ",
          "geometries.",
          call. = FALSE
        )
      }

      self$features <- features
    },

    # Recognise both terra (`SpatRaster`) and legacy raster package classes,
    # without requiring either package to be attached.
    is_raster = function(x) {
      inherits(x, c("SpatRaster", "RasterLayer", "RasterStack", "RasterBrick"))
    }
  )
)

#' @title Import a landscape image as a landscape object
#'
#' @description
#' Converts an image file into a [landscape] object: `sites` is stored as a
#' matrix of x/y pixel coordinates, and `features` is stored as a
#' tibble of greyscale intensity and/or RGBO colour channel values.
#'
#' @param file Path to an image file. Must be one of a png, jpeg, or bitmap file.
#' @param scale Scale the image up (>1) or down (<1).
#' @param keep_channels Logical. If `TRUE`, keep all colour channel data in
#'   `features` alongside the greyscale intensity. If `FALSE` (default),
#'   `features` contains only the greyscale intensity.
#'
#' @returns A [landscape] object with `sites` as an x/y coordinate matrix and
#'   `features` as a tibble of pixel intensity (and optionally colour
#'   channel) values.
#' @export
#'
#' @examples # See vignettes/generate-landscape-vignette.Rmd or eigenmove-vignette.Rmd for example usage
em_loadlandscape <- function(file,
                             scale = 1,
                             keep_channels = FALSE) {

  # load the image, resize it, and unfold it into one row per pixel with
  # one column per colour channel
  pixels <- imager::load.image(file) |>
    imager::imresize(scale = scale) |>
    as.data.frame(wide = "c") |>
    # flip so the image appears as it should when plotted
    dplyr::mutate(y = max(y) - y + 1) |>
    # greyscale intensity = mean across colour channels
    dplyr::rowwise() |>
    dplyr::mutate(intensity = mean(dplyr::c_across(dplyr::starts_with("c.")))) |>
    dplyr::ungroup()

  sites <- as.matrix(pixels[, c("x", "y")])
  colnames(sites) <- c("x", "y")

  if (keep_channels) {
    feature_cols <- c(grep("^c\\.", names(pixels), value = TRUE), "intensity")
  } else {
    feature_cols <- "intensity"
  }
  features <- dplyr::as_tibble(pixels[, feature_cols, drop = FALSE])

  landscape$new(sites = sites, features = features)
}
