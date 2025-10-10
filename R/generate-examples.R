## Functions related to loading and visualizing example landscapes for testing and demonstrations

#' @title Import an image file of an example landscape and convert it to a data frame
#'
#' @param file Path to an image file. Must be one of a png, jpeg, or bitmap file
#' @param scale Scale the image up (>1) or down (<1)
#' @param keep_channels Logical value to check if you want to keep all colour channel data, or just greyscale intensities
#'
#' @returns A dataframe specifying the coordinates (x and y), colour channel values, greyscale intensity, and habitat quality of each pixel
#' @export
#'
#' @examples # insert example
em_load_landscape <- function(file,
                             scale = 1,
                             keep_channels = FALSE
                             ){
  # create a new object by loading an external image file
  landscape = imager::load.image(file) |> # Use the path to an image of a landscape as the argument, e.g. png = "images/landscape.png"
    imager::imresize(scale = scale) |> # scale up (>1) or down (<1)
    # turn it into a data frame with one column per colour channel
    as.data.frame(wide = "c") |>
    #flips the image so that it appears as it should in the image file
    dplyr::mutate(y = max(y) - y + 1) |>
    # calculate a new intensity variable that represents the grey-scale value of the image
    dplyr::rowwise() |>
    dplyr::mutate(intensity = mean(dplyr::c_across(starts_with("c.")))) |>
    dplyr::ungroup() |>
    dplyr::mutate(habitat_type = dplyr::case_when(
      # assign qualities high mid low to pixels.
      intensity < 0.5 ~ "high",
      intensity > 0.9 ~ "low",
      TRUE ~ "mid"),
      habitat_type = factor(habitat_type, levels = c("low", "mid", "high")))

  if(!keep_channels){
    landscape <- dplyr::select(landscape,  !starts_with("c."))
  }
  return(landscape)
}


# Toy random walk movement model function.
# Calculates entries of the movement matrix for simple landscapes
calc_step <- function(dist, habitat_from, habitat_to,
                     step_length,
                     speed,
                     pref_strength){

  # These conditions end the function if something wonky is going on
  # stop if the distance _to_ is somehow different from the distance _from_
  stopifnot(length(dist) == length(habitat_from) &
              length(dist) == length(habitat_to))
  stopifnot(length(speed)==3 & all(speed>0))
  stopifnot(length(step_length)==3 & all(step_length>0))
  stopifnot(length(pref_strength)==3 & pref_strength>0)
  stopifnot(is.numeric(dist))   # stop if the distance is somehow not numeric
  stopifnot(is.character(habitat_to))   # stop if the habitat qualities are somehow numeric (or not characters)
  stopifnot(is.character(habitat_from))
  stopifnot(all(habitat_from %in% c("high", "mid", "low")))   # stop if the habitat qualities are anything but "high", "mid" or "low"
  stopifnot(all(habitat_to %in% c("high", "mid", "low")))
  from <- case_when(habitat_from =="high" ~ 1,
                    habitat_from =="mid" ~ 2,
                    habitat_from =="low" ~3)
  to <- case_when(habitat_to =="high" ~ 1,
                  habitat_to =="mid" ~ 2,
                  habitat_to =="low" ~ 3)
  # exponential function of Euclidean distance, scaled by habitat type of the leaving step
  base_step <- exp(-(dist-1)/step_length[from])

  # setting up the probability bandwidth parameter sigma,
  # higher probability of traveling to a higher quality habitat
  step_pref <- pref_strength[to]/pref_strength[from]
  step_speed <- speed[from]
  return(base_step*step_pref*step_speed)
}


# Generate a movement matrix given a landscape image input.
# Random walk model based on simple low, medium and high quality
# habitat distinction defined in image_to_dataframe
em_create_example_Q <- function(landscape,
                               step_length = c(0.5,0.5,2),
                               speed = c(0.5,0.5,2),
                               pref_strength = c(4,2,1)){ #argument is the landscape dataframe

  # Create an empty movement matrix with number of rows and columns each equal
  #to the total number of points on the landscape (number of rows in dataframe)
  n_pixels = nrow(landscape)
  movement_matrix = matrix(0, nrow = n_pixels, ncol = n_pixels)

  # Matrix whose entries are every pairwise combination of habitat qualities
  landscape_quality_matrix = outer(landscape$type,
                                   landscape$type,
                                   FUN = paste,
                                   sep = "-")

  # Matrix whose entries are every pairwise combination of Euclidean distances
  distance = as.matrix(dist(landscape[,c("x","y")]))

  # Core function to calculate entries of the movement matrix


  # For loop to run calc_step for each pair of points
  for(i in 1:n_pixels){
    # Ensure pairwise combinations of habitat quality types are in the correct
    # one-way order so that sigma properly takes into account the order of
    # "from" and "to" quality types
    quality_types = str_split_fixed(landscape_quality_matrix[,i],
                                    n = 2,
                                    pattern = "-")
    # Execute calc_step for each entry
    movement_matrix[,i] = calc_step(dist = distance[,i],
                                    habitat_from = quality_types[,2],
                                    habitat_to = quality_types[,1],
                                    step_length = step_length,
                                    speed = speed,
                                    pref_strength = pref_strength)

  }

  # Set up the diagonal of the movement_matrix to be something useful
  diag(movement_matrix) = -(colSums(movement_matrix)-diag(movement_matrix))
  return(movement_matrix)
}

# Simulating GP landscapes ####
# for use with simple random walk toy movement model
# 1 --- Run create_GP_landscape to generate a data frame.
#       e.g. current_GP <- create_GP_landscape()
# 2 --- Run rescale_landscape using the value column of the create_GP_landscape
#       to generate a new column with discrete low, mid, high values.
#       e.g. current_GP$type <- rescale_landscape(current_GP$value)
create_GP_landscape = function(landscape_width = 10,
                               landscape_height = 10,
                               patch_scale = 1){
  #This function randomly generates a new landscape using what's called a
  #Gaussian process; this assumes that the random value at each point in the
  #landscape follows a normal distribution, but that the random values for
  #points close to one another are correlated, so they have similar values.
  #Probably the best intro to GPs is here:
  #https://distill.pub/2019/visual-exploration-gaussian-processes/ although it
  #focuses more on using GPs for fitting data than for simulating data

  #landscape_width: specifies the width of the landscape in pixels
  #landscape_height: specifies the height of the landscape in pixels
  #patch_scale: specifies the size of the patch. Larger values of patch_scale
  #correspond to higher correlations between distant points, so larger (but less
  #common) patches

  #Checking arguments:
  if(length(landscape_width)>1 | landscape_width<0 | landscape_width%%1 !=0)
    stop("landscape width has to be a single positive integer")
  if(length(landscape_height)>1 | landscape_height<0 | landscape_height%%1 !=0)
    stop("landscape height has to be a single positive integer")
  if(length(patch_scale)>1 | patch_scale <0 )
    stop("patch_scale has to be a single positive number")

  #Creating landscape to output:
  n_patches = landscape_width*landscape_height
  landscape = crossing(x= 1:landscape_width,
                       y= 1:landscape_height)

  #Creates a distance matrix based on the landscape
  dist_mat = as.matrix(dist(landscape))

  #Generates the covariance matrix of the Gaussian process. This is a Matern
  #covariance function. The smoothness argument just results in somewhat
  #irregularly-shaped patches. The Matern function is from the fields package.
  cov_mat  = Matern(dist_mat, range=patch_scale, smoothness = 2.5)

  #simulates from the Gaussian process, using the mvrnorm function from the mgcv
  #package
  sim = mvrnorm(n=1,
                mu = rep(0, times=n_patches),
                Sigma = cov_mat)

  #adds that simulation to the landscape then returns the landscape to the user.
  landscape$intensity = as.vector(sim)

  return(landscape)
}


rescale_landscape = function(value,
                             good_hab_min = 1,
                             mid_hab_min  = 0.5
){
  #This function re-scales a continuous-valued landscape with a continuous set
  #of values to low, medium, and high values consistent with what we used for
  #the movement model, using the case_when function from the dplyr package.
  type = case_when(value>good_hab_min~"high",
                   value>mid_hab_min~"mid",
                   TRUE~"low")
  type = factor(type, levels = c("low", "mid", "high"))
  return(type)

}



# Other functions ####
# These functions need comments
em_neighbourdist <- function(locations, maxdist, nn = 100, ncores = 1){
  #uses the st_nn function to find the nn nearest neighbours of each point
  #that are within maxdist of it.

  maxdist <- maxdist
  nn_grid <- nngeo::st_nn(x = locations, y = locations,
                          sparse = TRUE,
                          maxdist = maxdist,
                          k = nn,
                          returnDist = TRUE,
                          parallel = ncores)
  n <- nrow(data)

  #Transforms the list returned in nn_grid into lists of indices of rows and
  #columns
  start <- list()
  end <- list()
  dists <- list()

  for(i in 1:n){
    n_vals <- length(nn_grid$nn[[i]])
    start[[i]] <-  rep(i, times=n_vals)
    end[[i]] <- nn_grid$nn[[i]]
    dists[[i]] <- nn_grid$dist[[i]]
  }

  start <- unlist(start)
  end <- unlist(end)
  dists <- unlist(dists)

  out <- Matrix::sparseMatrix(j = start,i = end, x = dists, giveCsparse = FALSE)

  if(any(Matrix::colSums(out)==0) | any(Matrix::rowSums(out)==0)) {
    warning("At least one location does not have any neighbours the given value of maxdist")
  }

  out
}

sparse_dispersemat <- function(nn_distmat,
                               patch_qual,
                               d0,
                               qual_bias,
                               dist_effect,
                               alpha,
                               lambda,
                               qual0,
                               dmax,
                               dmin = 1e-12){
  stopifnot(class(nn_distmat)[1]=="dgTMatrix")
  n <- nrow(nn_distmat)
  n_nonzero <- length(nn_distmat@i)
  stopifnot(length(patch_qual)==n)

  disp_mat <- Matrix::sparseMatrix(i = nn_distmat@i+1,
                                   j = nn_distmat@j+1,
                                   x = 1,
                                   dims = c(n,n))

  #have to add 1 to indices as the dgTmatrix format starts indices at 0
  i_vals <- nn_distmat@i+1
  j_vals <- nn_distmat@j + 1
  dists <- nn_distmat@x

  start_qual <- patch_qual[j_vals]
  end_qual <- patch_qual[i_vals]

  base_rate <- d0 + d0*lambda*(plogis(-(start_qual-qual0)*alpha))
  val <-  base_rate*exp(qual_bias*(end_qual-start_qual))*exp(-dist_effect*dists)
  val <- ifelse(val>dmax, dmax, val)
  #always some tiny, but non-zero movement to all connected locations
  val <- ifelse(val<dmin, dmin, val)
  val[i_vals==j_vals] <- 0

  disp_mat <- Matrix::sparseMatrix(i = i_vals,
                                   j = j_vals,
                                   x = val,
                                   dims = c(n,n))

  diag(disp_mat) <- - (Matrix::colSums(disp_mat) - diag(disp_mat))
  disp_mat
}
