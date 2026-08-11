# Functions related to loading and visualizing example landscapes for testing and demonstrations

#' @title Generate a movement matrix given a landscape image input and step parameters
#'
#' @description
#' Generates a movement matrix based on a simple random walk model where step length, speed, and preference strength can be adjusted according to habitat quality types defined in the landscape data frame.
#'
#' @param landscape data frame of landscape pixel coordinates and habitat quality types
#' @param step_length vector of step lengths for each habitat quality type
#' @param speed vector of step speeds for each habitat quality type
#' @param pref_strength vector of preference strengths for each habitat quality type
#'
#' @returns A symmetrical movement matrix: A square matrix giving the probability of movement from any point of the landscape to any other point. Size is the number of points in the landscape.
#' @export
#'
#' @examples See vignettes/generate-landscape-vignette.Rmd for example usage
em_create_example_Q <- function(landscape,
                               step_length = c(0.5,0.5,2),
                               speed = c(0.5,0.5,2),
                               pref_strength = c(4,2,1)){ #argument is the landscape data frame

  # Create an empty movement matrix with number of rows and columns each equal
  #to the total number of points on the landscape (number of rows in data frame)
  n_pixels = nrow(landscape)
  movement_matrix = matrix(0, nrow = n_pixels, ncol = n_pixels)

  # Matrix whose entries are every pairwise combination of habitat qualities
  landscape_quality_matrix = outer(landscape$type,
                                   landscape$type,
                                   FUN = paste,
                                   sep = "-")

  # Matrix whose entries are every pairwise combination of Euclidean distances
  distance = as.matrix(stats::dist(landscape[,c("x","y")]))

  # Core function to calculate entries of the movement matrix


  # For loop to run calc_step for each pair of points
  for(i in 1:n_pixels){
    # Ensure pairwise combinations of habitat quality types are in the correct
    # one-way order so that sigma properly takes into account the order of
    # "from" and "to" quality types
    quality_types = stringr::str_split_fixed(landscape_quality_matrix[,i],
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


#' @title Simulate a toy landscape using a Gaussian process
#'
#' @description
#' Simulate GP landscapes for use with simple random walk toy movement models. A Gaussian process  assumes that a random value (habitat quality) at each point in a landscape follows a normal distribution, and that the random values for points close to one another are correlated, so they have similar values.
#'
#' @param landscape_width width of image in pixels
#' @param landscape_height height of image in pixels
#' @param patch_scale size of the patch. Larger values of patch_scale correspond to higher correlations between distant points, so larger (and fewer) patches
#'
#' @returns A data frame specifying the coordinates (x and y) and intensity value of each pixel
#' @export
#'
#' @examples See vignettes/generate-landscape-vignette.Rmd for example usage
create_GP_landscape = function(landscape_width = 10,
                               landscape_height = 10,
                               patch_scale = 1){

  #Checking arguments:
  if(length(landscape_width)>1 | landscape_width<0 | landscape_width%%1 !=0)
    stop("landscape width has to be a single positive integer")
  if(length(landscape_height)>1 | landscape_height<0 | landscape_height%%1 !=0)
    stop("landscape height has to be a single positive integer")
  if(length(patch_scale)>1 | patch_scale <0 )
    stop("patch_scale has to be a single positive number")

  #Creating landscape to output:
  n_patches = landscape_width*landscape_height
  landscape = tidyr::crossing(x= 1:landscape_width,
                              y= 1:landscape_height)

  #Creates a distance matrix based on the landscape
  dist_mat = as.matrix(stats::dist(landscape))

  #Generates the covariance matrix of the Gaussian process. This is a Matern
  #covariance function. The smoothness argument just results in somewhat
  #irregularly-shaped patches. The Matern function is from the fields package.
  cov_mat  = fields::Matern(dist_mat, range=patch_scale, smoothness = 2.5)

  #simulates from the Gaussian process, using the mvrnorm function from the mgcv
  #package
  sim = MASS::mvrnorm(n=1,
                mu = rep(0, times=n_patches),
                Sigma = cov_mat)

  #adds that simulation to the landscape then returns the landscape to the user.
  landscape$intensity = as.vector(sim)

  return(landscape)
}


#' @title Rescale continuous-valued landscape to discrete habitat quality types
#'
#' @description
#' This function re-scales a continuous-valued landscape with a continuous set of values to low, medium, and high values consistent with what is used for the movement model (`em_create_example_Q()`), using the case_when function from the dplyr package.
#'
#' @param value Continuous-valued landscape variable
#' @param good_hab_min Threshold value for high quality habitat
#' @param mid_hab_min Threshold value for medium quality habitat
#'
#' @returns A factor vector with levels "low", "mid", and "high"
#' @export
#'
#' @examples See vignettes/generate-landscape-vignette.Rmd for example usage
rescale_landscape = function(value,
                             good_hab_min = 1,
                             mid_hab_min  = 0.5){

  type = dplyr::case_when(value>good_hab_min~"high",
                   value>mid_hab_min~"mid",
                   TRUE~"low")

  type = factor(type, levels = c("low", "mid", "high"))

  return(type)
}


#' @title Generate a sparse matrix of nearest neighbour distances
#'
#' @description
#' This function uses the st_nn function from the nngeo package to find the nearest neighbours of each point within a maximum distance, and returns a sparse matrix of the distances to those neighbours.
#'
#' @param locations An sf or sfc object with point geometries
#' @param maxdist Maximum distance to consider for neighbours
#' @param nn Number of nearest neighbours to consider
#' @param ncores Number of cores to use for parallel processing
#'
#' @returns A sparse matrix of nearest neighbour distances
#' @export
#'
#' @examples See vignettes/generate-landscape-vignette.Rmd for example usage
em_neighbourdist <- function(locations, maxdist, nn = 100, ncores = 1){

  maxdist <- maxdist
  n <- nrow(utils::data)
  nn_grid <- nngeo::st_nn(
    x = locations, y = locations, # sf coords
    sparse = TRUE,
    maxdist = maxdist,
    k = nn, # set to n (number of landscape points) to get all within maxdist
    returnDist = TRUE,
    parallel = ncores)

  #Transforms the list returned in nn_grid into lists of indices of rows and columns
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

# This function needs more comments and better descriptions of parameters

#' @title Generate a sparse dispersal matrix
#'
#' @description
#' This function generates a sparse dispersal matrix based on a nearest neighbour distance matrix and patch qualities.
#'
#' @param nn_distmat Sparse matrix of nearest neighbour distances
#' @param patch_qual Vector of patch qualities
#' @param d0 Base dispersal rate
#' @param qual_bias Quality bias parameter
#' @param dist_effect Distance effect parameter
#' @param alpha Parameter used in calculation of base movement rate
#' @param lambda Parameter used in calculation of base movement rate
#' @param qual0 Parameter used in calculation of base movement rate
#' @param dmax Maximum movement rate
#' @param dmin Minimum movement rate
#'
#' @returns A sparse dispersal matrix
#' @export
#'
#' @examples TBD
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

  base_rate <- d0 + d0*lambda*(stats::plogis(-(start_qual-qual0)*alpha))
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
