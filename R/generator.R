#' Generator Class
#'
#' @description
#' Continuous-time Markov chain describing rates of transition between points on a landscape.
#'
#' @export
Generator <- R6::R6Class(
  "Generator",
  public = list(
    #' @field id Generator identifier
    id = NULL,

    #' @field landscape Reference to a Landscape object
    landscape = NULL,

    #' @description Create a new generator
    #' @param id Generator identifier
    #' @param landscape A Landscape object
    initialize = function(id, landscape) {
      self$id <- id
      self$landscape <- landscape
    }
  )
)
