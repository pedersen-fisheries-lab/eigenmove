# New Eigenmove Functions ####

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

#' @title Build a random walk movement model
#'
#' @description
#' Constructs a movement model function based on a simple random walk model.
#' Step length, step speed, and habitat preference strength can be adjusted
#' according to habitat quality type. This functional returns a `movemodel`
#' function that can later be applied to any landscape (e.g. via
#' `em_buildgenerator()`) to produce a generator matrix.
#'
#' @param step_length vector of step lengths for each habitat quality type
#' @param step_speed vector of step speeds for each habitat quality type
#' @param pref_strength vector of preference strengths for each habitat quality type
#'
#' @returns A function (`movemodel`) that,when applied to a `landscape` via
#' `em_buildgenerator` returns a generator matrix describing movement rates on
#' that `landscape`
#' @export
#'
#' @examples
#' movemodel <- em_randomwalk()
em_randomwalk <- function(step_length = c(0.5, 1, 2),
                          step_speed = c(0.5, 1, 2),
                          pref_strength = c(4, 2, 1)) {

  stopifnot(length(step_length) == 3 & all(step_length > 0))
  stopifnot(length(step_speed) == 3 & all(step_speed > 0))
  stopifnot(length(pref_strength) == 3 & all(pref_strength > 0))

  movemodel <- function(landscape) {

    n_pixels <- nrow(landscape)
    dispersal_matrix <- matrix(0, nrow = n_pixels, ncol = n_pixels)

    landscape_quality_matrix <- outer(landscape$type,
                                      landscape$type,
                                      FUN = paste,
                                      sep = "-")

    distance <- as.matrix(stats::dist(landscape[, c("x", "y")]))

    for (i in 1:n_pixels) {
      quality_types <- stringr::str_split_fixed(landscape_quality_matrix[, i],
                                                n = 2,
                                                pattern = "-")
      dispersal_matrix[, i] <- calc_step(dist = distance[, i],
                                        habitat_from = quality_types[, 2],
                                        habitat_to = quality_types[, 1],
                                        step_length = step_length,
                                        step_speed = step_speed,
                                        pref_strength = pref_strength)
    }

    diag(dispersal_matrix) <- -(colSums(dispersal_matrix) - diag(dispersal_matrix))
    dispersal_matrix
  }
  class(movemodel) <- c("movemodel", class(movemodel))
  return(movemodel)
}

#' @title Toy random walk step function
#'
#' @description
#' Performs random walk steps on landscape according to `movemodel`
#'
#'
#' @param dist Distance matrix. A matrix whose entries are every pairwise combination of Euclidean distances
#' @param habitat_from type of habitat the step is starting from
#' @param habitat_to type of habitat the step is going to
#' @param step_length size of step
#' @param step_speed speed of step
#' @param pref_strength preference strength for habitat types
#'
#' @returns A vector of movement probabilities
#' @export
#'
#' @examples Used internally by `movemodel`
calc_step <- function(dist,
                      habitat_from,
                      habitat_to,
                      step_length,
                      step_speed,
                      pref_strength){

  # These conditions end the function if something wonky is going on
  # stop if the distance _to_ is somehow different from the distance _from_
  stopifnot(length(dist) == length(habitat_from) &
              length(dist) == length(habitat_to))
  stopifnot(length(step_speed)==3 & all(step_speed>0))
  stopifnot(length(step_length)==3 & all(step_length>0))
  stopifnot(length(pref_strength)==3 & pref_strength>0)
  stopifnot(is.numeric(dist))   # stop if the distance is somehow not numeric
  stopifnot(is.character(habitat_to))   # stop if the habitat qualities are somehow numeric (or not characters)
  stopifnot(is.character(habitat_from))
  stopifnot(all(habitat_from %in% c("high", "mid", "low")))   # stop if the habitat qualities are anything but "high", "mid" or "low"
  stopifnot(all(habitat_to %in% c("high", "mid", "low")))

  # Weights
  from <- dplyr::case_when(habitat_from =="high" ~ 1,
                           habitat_from =="mid" ~ 2,
                           habitat_from =="low" ~3)
  to <- dplyr::case_when(habitat_to =="high" ~ 1,
                         habitat_to =="mid" ~ 2,
                         habitat_to =="low" ~ 3)

  # exponential function of Euclidean distance, scaled by habitat type of the leaving step
  base_step <- exp(-(dist-1)/step_length[from])

  # setting up the probability bandwidth parameter sigma,
  # higher probability of traveling to a higher quality habitat
  step_pref <- pref_strength[to]/pref_strength[from]
  step_speed <- step_speed[from]
  return(base_step*step_pref*step_speed)
}

#' @title Build a generator matrix from a movement model and a landscape
#'
#' @description
#' Applies a `movemodel` function (e.g. produced by `em_randomwalk()`) to a
#' specific landscape to produce a generator matrix.
#'
#' @param landscape data frame of landscape pixel coordinates and habitat quality types
#' @param movemodel a movement model function that describes individual
#' movement on a landscape
#'
#' @returns A symmetrical CTMC generator matrix: a square matrix giving the
#' probability (or rate) of movement from any point of the landscape to any
#' other point. Size is the number of points in the landscape.
#' @export
#'
#' @examples
#' movemodel <- em_randomwalk()
#' # generator <- em_buildgenerator(landscape, movemodel)
em_buildgenerator <- function(landscape, movemodel) {
  stopifnot(is.function(movemodel))
  generator <- movemodel(landscape)
  class(generator) <- c("generator", class(generator))

  structure(list(landscape = landscape, generator = generator),
            class = "generator")
}

em_loadgenerator <- function(landscape, generator){

}

em_simmove <- function(generator,
                       steps = 5,
                       replicates = 3,
                       origins = 2) {

  landscape <- generator$landscape
  generator <- generator$generator
  n_locs <- nrow(landscape)

  # draw `origins` random starting locations from the landscape
  origin_samples <- sample(seq_len(n_locs), size = origins, replace = TRUE)
  paths <- vector("list", origins)
  names(paths) <- paste0("origin_", seq_len(origins))

  # function to simulate one Markov-chain path of `steps` steps from origin
  sample_path <- function(origin, steps, generator) {
    path <- integer(steps + 1) # each path is origin and `steps` steps
    path[1] <- origin
    for (s in seq_len(steps)) { # take steps according to probabilities in generator
      probs <- generator[, path[s]] # column arranged stochastic matrix
      probs[path[s]] <- 0 # don't include "stay in place"
      path[s + 1] <- sample(seq_along(probs), size = 1, prob = probs)
    }
    path
  }

  # main loop. Repeat `origins` number of times
  for (i in seq_len(origins)) {
    # Get the origin location generated earlier
    origin_xy <- landscape[origin_samples[i], c("x", "y")]

    # Set up replicates list to nest within paths
    rep_list <- vector("list", replicates)
    names(rep_list) <- paste0("rep_", seq_len(replicates))

    # Generate a replicate for path/origin i
    for (r in seq_len(replicates)) {
      loc_seq <- sample_path(origin_samples[i], steps, generator) # take steps
      rep_list[[r]] <- data.frame( # organize replicates list
        step      = 0:steps,
        loc_index = loc_seq,
        x         = landscape$x[loc_seq],
        y         = landscape$y[loc_seq]
      )
    }
    # Store replicates nested within each origin
    paths[[i]] <- list(
      origin       = c(x = origin_xy$x, y = origin_xy$y),
      origin_index = origin_samples[i],
      replicates   = rep_list
    )
  }
  paths
}

em_probmove <- function(x, ...) {
  UseMethod("em_probmove")
}

em_probmove.generator <- function(x, ...) {
  # x is your generator object
  cat("Processing generator object\n")

  # Your specific implementation for generator class
  # For example:
  result <- list(
    type = "generator",
    data = x,
    steps = steps,
    timestamp = Sys.time()
  )
  class(result) <- "em_probmove_result"
  return(result)
}

em_probmove.eigenmove <- function(x, replicates = 2, ...) {
  cat("Processing eigenmove object\n")

  # Your specific implementation for eigenmove class
  result <- list(
    type = "eigenmove",
    data = x,
    replicates = replicates,
    timestamp = Sys.time()
  )
  class(result) <- "em_probmove_result"
  return(result)
}
