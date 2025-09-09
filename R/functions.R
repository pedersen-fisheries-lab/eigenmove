# Compiling c++ distance functions ####

# Eigenpatch Functions  ####

# Import and set up landscape as a data frame usable by calc_step.
#' Title
#'
#' @param png An image file
#' @param scale Scale the image up (>1) or down (<1)
#'
#' @returns A dataframe specifying the habitat quality of each pixel
#' @export
#'
#' @examples # insert example
image_to_dataframe = function(png, scale=1){
  # create a new object by loading an external image file (%>% = and then...)

  landscape = load.image(png) %>% # Use the path to an image of a landscape as the argument, e.g. png = "images/landscape.png"
    imresize(scale = scale) %>% # scale up (>1) or down (<1)
    as.data.frame() %>%    # turn it into a dataframe
    group_by(x,y) %>%    # temporarily group each pixel
    summarize(value = mean(value), .groups = 'drop') %>%    # make each pixel the mean of its RGBO values, then ungroup
    mutate(type = case_when(     # assign qualities high mid low to pixels.
      value < 0.5 ~ "high",
      value ==  1 ~ "low",
      TRUE ~ "mid"),
      type = factor(type, levels = c("low", "mid", "high")))
  return(landscape)
}

# Toy random walk movement model function.
# Calculates entries of the movement matrix for simple landscapes
calc_step = function(dist, habitat_from, habitat_to,
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


# Generate a dispersal matrix given a landscape image input.
# Random walk model based on simple low, medium and high quality
# habitat distinction defined in image_to_dataframe
calculate_dispersal = function(landscape,
                               step_length = c(0.5,0.5,2),
                               speed = c(0.5,0.5,2),
                               pref_strength = c(4,2,1)){ #argument is the landscape dataframe

  # Create an empty dispersal matrix with number of rows and columns each equal
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

  # Core function to calculate entries of the dispersal matrix


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

check_conjugate_state = function(eigenvalues,tol = 1e-14){
  #this function is used to check that, if the smallest eigenvalue has a
  #non-zero imaginary component, and if so, if that value is matched with a
  #conjugate eigenvalue in the next smallest value. Otherwise, the kinetic
  #distance calculation won't result in all real kinetic distances
  stopifnot(is.numeric(eigenvalues)|is.complex(eigenvalues))
  if(abs(Im(eigenvalues[1]))<tol){
    val = TRUE
  } else if(dplyr::near(eigenvalues[1]*eigenvalues[2] ,
                        abs(eigenvalues[1])^2,
                        tol = tol)){
      val = TRUE
  } else{
    val = FALSE
  }
  val
}

# Calculate d+1 eigenvalues and left and right eigenvectors.
calculate_eigenfunctions = function(movement_matrix, d,  sigma_val = 1e-8){
  #d is the number of dimensions to retain if sigma_val is a non-null number,
  #uses the "shift-and-invert" mode to calculate eigenvectors and values. It
  #works better for very large matrices. Set to NULL if you don't want to use
  #this method


  right_eigs = eigs(movement_matrix, k= d, which = "LM",sigma = sigma_val)
  left_eigs = eigs(t(movement_matrix), k= d, which = "LM",sigma = sigma_val)


  stopifnot(all.equal(abs(right_eigs$values[-d]/left_eigs$values[-d]),
                      rep(1, times = d-1),
                      tolerance = 1e-7))


  #ensure that the left eigenvectors are scaled appropriately so the left
  #eigenvectors form the inverse of the right eigenvector matrix
  for(i in 1:d){
    q_val = sum(right_eigs$vectors[,i]*left_eigs$vectors[,i])
    left_eigs$vectors[,i] = left_eigs$vectors[,i]/q_val
  }


  eigenfunctions_list = list(Phi = right_eigs$vectors,
                             Psi = left_eigs$vectors,
                             Lambda = right_eigs$values,
                             d = d)
  return(eigenfunctions_list)
}

# Calculate a matrix of diffusion distances based on the movement matrix and eigenfunctions.
calculate_kinetic_distances = function(movement_matrix,
                                       d,
                                       T,
                                       sigma_val=1e-8,
                                       discrete_time = FALSE,
                                       scale_by_density = FALSE,
                                       keep_imaginary = FALSE,
                                       progress_bar = FALSE){
  stopifnot(length(T)==1) # T needs to be one number, not a vector
  stopifnot(T>0) # T needs to be positive
  stopifnot(length(d)==1) # number of dimensions specified needs to be one number, not a vector
  stopifnot(d>1) # number of dimensions needs to be more than 1
  if(d%%1 != 0) stop("d must be a positive integer greater than 1") # checks remainder when 1 divides d
  if(discrete_time &  T%%1 != 0){
    stop("If using a discrete-time dispersal matrix, the time scale T must be an integer")
  }
  if(discrete_time & !all(near(colSums(movement_matrix),y = 1,tol = 1e-10))){
    stop("If using a discrete-time random walk, the columns of the movement matrix must sum to one, and all entries must be positive")
  } else if(!all(near(colSums(movement_matrix),y=0, tol=1e-10))){
    stop("If using a continuous-time random walk, the columns of the movement matrix must sum to zero")
  }

  n_pixels = nrow(movement_matrix) # number of pixels

  # Left and right eigenfunctions, need d+1 because the leading right eigenvector
  # is long term distribution and considered separately for diffusion distances
  eigendecomp = calculate_eigenfunctions(movement_matrix = movement_matrix,
                                         d = d+1,
                                         sigma_val = sigma_val)

  #convert the eigenvalues of the eigenvector decomposition to their exponential
  #values, scaling by T.
  if(discrete_time){
    if((1- abs(eigendecomp$Lambda[d])) < 1e-15){
      stop("Landscape is not fully connected (more than one eigenvalue of the movement matrix is equal to 1)")
    }
    timescale_matrix <- outer(eigendecomp$Lambda[1:d], eigendecomp$Lambda[1:d],
                              FUN = "*")
    timescale_matrix <- (timescale_matrix - timescale_matrix^(T + 1))/
      (1 - timescale_matrix)
  }else{
    # If there's more than one zero eigenvalue, the patches aren't connected,
    # stop the program and alert the user
    if(abs(eigendecomp$Lambda[d]) < 1e-15){
      stop("Landscape is not fully connected (more than one eigenvalue of the movement matrix is equal to 0)")
    }
    #this is based off taking the average exponentiated value across time length
    #of time T
    timescale_matrix <- outer(eigendecomp$Lambda[1:d], eigendecomp$Lambda[1:d],
                              FUN = "+")
    timescale_matrix <- (exp(T*timescale_matrix)-1)/(timescale_matrix*T)
  }

  # If (when exponentiated and scaled by T) the eigenvalue with the smallest
  # real part (1st lambda) is more than 5% of the eigenvalue with the largest
  # real part (dth lambda), alert the user. More eigenvectors/eigenvalues are
  # needed to capture fine-grained movement detail.
  if(discrete_time & Re(exp(eigendecomp$Lambda[d]) / exp(eigendecomp$Lambda[1])) > 0.05) {
    warning("n_eigs potentially too low to capture fine-grained movement details. Consider increasing n_eigs")
  }


  occupancy_prob = Re(eigendecomp$Phi[,d+1]) # Stable distribution of landscape, leading right eigenvector
  occupancy_prob = occupancy_prob/sum(occupancy_prob) # Standardized

  eigendecomp$Lambda = eigendecomp$Lambda[-(d+1)] # Vector of Eigenvalues
  eigendecomp$Phi = eigendecomp$Phi[,-(d+1)] # Matrix of right eigenvectors
  eigendecomp$Psi = eigendecomp$Psi[,-(d+1)] # Matrix of left eigenvectors

  if(!check_conjugate_state(eigendecomp$Lambda)){
    #check to see if the last eigenvalue has a matching complex conjugate (if complex)
    #if not, drop the last eigenvalue
    if(!check_conjugate_state(eigendecomp$Lambda[-1])){
      #down the line, need to figure out a fix for when Rspectra occasionally
      #does not return complex conjugates. For now I'll leave this warning in
      #place
      stop("At least one of the eigenvalues does not have a complex conjugate")
    }
    d = d-1
    eigendecomp$d = d
    eigendecomp$Psi  =  eigendecomp$Psi[,-1]
    eigendecomp$Phi  =  eigendecomp$Phi[,-1]
    eigendecomp$Lambda = eigendecomp$Lambda[-1]

    warning(paste0("Decreased the number of eigenvectors used by 1 to ", d, " as the final eigenvalue did not have a complex conjugate for the specified d value"))
  }

  # If using DBSCAN for clustering, this gives a lower limit for
  # the epsilon parameter.
  eigendecomp$eps_threshold <- abs(
    ((eigendecomp$Lambda %*% eigendecomp$Lambda) -
       (eigendecomp$Lambda^(T+1) %*% eigendecomp$Lambda^(T+1))) /
      (1 - (eigendecomp$Lambda %*% eigendecomp$Lambda)))

  # Calculates a matrix of inner products of the right eigenvectors either
  # scaled or not scaled by the inverse of patch-specific long-term occupancy
  if(scale_by_density){
    inv_occupancy <- diag(1/occupancy_prob)
    inner_matrix <- t(eigendecomp$Phi)%*%inv_occupancy^2 %*%eigendecomp$Phi

  } else{
    inner_matrix <-  t(eigendecomp$Phi) %*% eigendecomp$Phi
  }

  #elementwise (Kroenecker-product) of the Phi and timescale matrices
  inner_matrix <- inner_matrix*timescale_matrix

  diff_list <- list()
  dists <- fastdist(eigendecomp$Psi, inner_matrix)
  gc()

  #That calculated the squared diffusion distances; we return the unsquared values
  if(!keep_imaginary){
    dists = Re(dists)
  }
  dists = sqrt(dists)

  attributes(dists) <- list(method = "kinetic", # Give output object the class 'dist' (distance matrix)
                           Diag = FALSE,
                           Upper = FALSE,
                           Size = n_pixels,
                           class = "dist")
  out <- list(dists = dists, # diffusion distances as dist object
              occupancy_prob = occupancy_prob, # long term distribution
              eps_threshold = eigendecomp$eps_threshold,
              d = d,
              T = T)
  return(out)
}

calculate_clusters = function(cluster_type = c("hclust", "DBSCAN", "OPTICS"),
                              landscape,
                              out,
                              min_dens = 1/nrow(landscape),
                              n_clust = n_clust,
                              ...){
  parms <- list(...)
  cluster_type = match.arg(cluster_type) # associate w/ argument " "

  # Set-up for all clustering algorithms
  landscape$clusters = NA # init empty clusters column, going to fill only in_patch entries
  landscape$dens = out$occupancy_prob # long term occupancy density
  landscape$in_patch = landscape$dens > min_dens # definition of in_patch points. Modify argument min_dens in call to function to change threshold
  cluster_setup = as.matrix(out$dists)[landscape$in_patch,landscape$in_patch] # need to read as matrix to subset in_patch points for clustering
  cluster_setup = as.dist(cluster_setup) # back to a distance object

  cluster_type = match.arg(cluster_type) # associate w/ argument " "

  landscape = landscape %>%   # Add columns for clusters, density, and in_patch. min_dens is threshold
    mutate(clusters = NA,
           dens = out$occupancy_prob,
           in_patch = dens > min_dens)

  if(cluster_type == "hclust") { # Hierarchical agglomerative clustering
    hclust_clusters <- hclust(cluster_setup, ...)
    hclust_clusters <- cutree(hclust_clusters, k = n_clust)
    landscape$clusters[landscape$in_patch] <- hclust_clusters
  }

    if(cluster_type == "DBSCAN") {
    if(eps < out$eps_threshold) {
      stop(paste0("Based on landscape size and structure, eps must be larger than ",
                  round(out$eps_threshold, digits = 2),
                  ". See appendix for justification. For general DBSCAN parameter setting guidelines see DBSCAN documentation"))
    }
    dbscan_clusters <- dbscan(cluster_setup, ...)$cluster
    landscape$clusters[landscape$in_patch] <- dbscan_clusters
  }

  if(cluster_type == "OPTICS") {
    optics_clusters <- do.call(optics, list(x = cluster_setup, minPts = parms$minPts))
    reachability <- optics_clusters
    optics_clusters <- do.call(extractDBSCAN, list(object = optics_clusters, eps_cl = parms$eps_cl))
    optics_clusters <- optics_clusters$cluster
    landscape$clusters[landscape$in_patch] <- optics_clusters
    optics_out <- list("landscape" = landscape, "reachability" = reachability)
    return(optics_out)
    break
  }
  return(landscape)
}

# Simulating GP landscapes ####
# for use with simple random walk toy movement model
# 1 --- Run create_GP_landscape to generate a data frame.
#       e.g. current_GP <- create_GP_landscape()
# 2 --- Run rescale_landscape using the value column of the create_GP_landscape
#       to generate a new column with discrete low, mid, high values.
#       e.g. current_GP$type <- rescale_landscape(current_GP$value)
create_GP_landscape = function(landscape_width=10,
                                  landscape_height=10,
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
  landscape$value = as.vector(sim)

  return(landscape)
}


rescale_landscape = function(value,
                             good_hab_min = 1,
                             mid_hab_min  = 0.5
                             ){
  #This function re-scales a continuous-valued landscape with a continuous set
  #of values to low, medium, and high values consistent with what we used for
  #the dispersal model, using the case_when function from the dplyr package.
  type = case_when(value>good_hab_min~"high",
                        value>mid_hab_min~"mid",
                        TRUE~"low")
  type = factor(type, levels = c("low", "mid", "high"))
  return(type)

}

# Other functions ####
# These functions need comments
sparse_distmat <- function(data, maxdist, nn = 100, ncores=1){
  #uses the st_nn function to find the nn nearest neighbours of each point
  #that are within maxdist of it.

  maxdist <- maxdist
  nn_grid <- nngeo::st_nn(data,data,
                          sparse=TRUE,
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

  out <- Matrix::sparseMatrix(j = start,i =end, x = dists,giveCsparse = FALSE)
  if(any(Matrix::colSums(out)==0) | any(Matrix::rowSums(out)==0)) warning("At least one location is disconnected at this maxdist")
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
  #always some tiny, but non-zero dispersal to all connected locations
  val <- ifelse(val<dmin, dmin, val)
  val[i_vals==j_vals] <- 0

  disp_mat <- Matrix::sparseMatrix(i = i_vals,
                                   j = j_vals,
                                   x = val,
                                   dims = c(n,n))


  diag(disp_mat) <- - (Matrix::colSums(disp_mat) - diag(disp_mat))
  disp_mat
}

calc_density <- function(disperse_mat, sigma=1e-16){
  dens = RSpectra::eigs(disperse_mat,k = 1,which = "LM",sigma = 1e-16)
  dens = Re(dens$vectors[,1])
  dens <- dens/sum(dens)
  dens
}

optimize_dispersal_model <- function(dist_mat,
                                     quality_index,
                                     ref_index,
                                     dist_effect,
                                     d0,
                                     dmax,
                                     alpha_start,
                                     qual0_start,
                                     bias_start,
                                     lambda_start){
  optim_func <- function(p){
    alpha <- p[1]
    qual0 <- p[2]
    bias <- p[3]
    lambda <- p[4]
    cat(p)
    disp_matrix <- sparse_dispersemat(nn_distmat = dist_mat,
                                      patch_qual = quality_index,
                                      d0 = d0,
                                      dist_effect = dist_effect,
                                      qual_bias = bias,
                                      alpha = alpha,
                                      lambda = lambda,
                                      qual0 = qual0,
                                      dmax = dmax)
    if(any(is.na(disp_matrix))) {
      value <- 1
    } else{
      current_density <- calc_density(disp_matrix)
      value <- 1-cor(current_density, ref_index)
    }
    cat(" ",value,"\n")
    gc()
    value
  }
  p_start <- c(alpha = alpha_start,
               qual0=qual0_start,
               bias=bias_start,
               lambda=lambda_start)
  optim_vals <- nlm(optim_func,p = p_start,stepmax = 1,ndigit = 3, fscale = 0)
  optim_vals

}
