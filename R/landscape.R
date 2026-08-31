#' @title Import a landscape as a data frame
#'
#' @description
#' Converts image file into a data frame with columns for x and y pixel cartesian coordinates, greyscale intensities and optionally RGBO color channels
#'
#' @param file Path to an image file. Must be one of a png, jpeg, or bitmap file
#' @param scale Scale the image up (>1) or down (<1)
#' @param keep_channels Logical value to check if you want to keep all colour channel data, or just greyscale intensities
#'
#' @returns A data frame specifying the coordinates (x and y) and pixel values
#' @export
#'
#' @examples # See vignettes/generate-landscape-vignette.Rmd or eigenmove-vignette.Rmd for example usage
em_loadlandscape <- function(file,
                             scale = 1,
                             keep_channels = FALSE){
  # create a new object by loading an external image file
  landscape = imager::load.image(file) |> # Use the path to an image of a landscape as the argument, e.g. png = "images/landscape.png"
    imager::imresize(scale = scale) |> # scale up (>1) or down (<1)
    # turn it into a data frame with one column per colour channel
    as.data.frame(wide = "c") |>
    #flips the image so that it appears as it should in the image file
    dplyr::mutate(y = max(y) - y + 1) |>
    # calculate a new intensity variable that represents the grey-scale value of the image
    dplyr::rowwise() |>
    dplyr::mutate(intensity = mean(dplyr::c_across(starts_with("c."))))

  if(!keep_channels){
    landscape <- dplyr::select(landscape,  !starts_with("c."))
  }
  class(landscape) <- c("landscape", class(landscape))
  return(landscape)
}

#' Landscape Class
#'
#' @description
#' Object type representing a landscape image in eigenmove
#'
#' @export
Landscape <- R6::R6Class(
  "Landscape",
  public = list(
    #' @field data Data frame with landscape data
    data = NULL,

    #' @field name Landscape name/identifier
    name = NULL,

    #' @description Create a new landscape object
    #' @param data Data frame with landscape data (must have x, y, intensity columns)
    #' @param name Optional name for the landscape
    initialize = function(data, name = NULL) {
      self$data <- data
      self$name <- name %||% "unnamed_landscape"
    },
  )
)
