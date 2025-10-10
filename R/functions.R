# Eigenmove Functions  ####



#' Check conjugate state of eigenvalues
#'
#'  This function checks that if the smallest eigenvalue has a
#'  non-zero imaginary component, that imaginary component is matched with a
#'  conjugate eigenvalue in the next smallest value. Otherwise, the kinetic
#'  distance calculation won't result in all real kinetic distances.
#'
#' @param eigenvalues A vector of eigenvalues
#' @param tol tolerance defining which numbers are considered zero
#'
#' @returns TBD
#' @export
#'
#' @examples TBD
check_conjugate_state = function(eigenvalues,tol = 1e-14){
  # This function checks that if the smallest eigenvalue has a
  #   non-zero imaginary component, that imaginary component is matched with a
  #  conjugate eigenvalue in the next smallest value. Otherwise, the kinetic
  #  distance calculation won't result in all real kinetic distances.

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

#' Title Calculate eigenvalues and eigenvectors of a movement matrix
#'
#' Calculate d+1 eigenvalues and left and right eigenvectors.
#'
#' @param movement_matrix Matrix giving the probability of movement
#' from any point of the landscape to any other point. A square matrix
#' whose size is the number of points in the landscape.
#' @param d Number of dimensions to retain
#' @param sigma_val If sigma_val is non-null, this function uses the "shift-and-invert"
#' mode to calculate eigenvectors and values. Recommended for very large matrices.
#' Set to NULL if you don't want to use this method
#'
#' @returns TBD
#' @export
#'
#' @examples TBD
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

#' Title Calculate kinetic distances
#'
#' Calculate a matrix of kinetic distances based on the movement matrix and eigenfunctions.
#'
#' @param movement_matrix Matrix giving the probability of movement
#' from any point of the landscape to any other point. A square matrix
#' whose size is the number of points in the landscape.
#' @param d Number of dimensions to retain
#' @param sigma_val If sigma_val is non-null, this function uses the "shift-and-invert"
#' mode to calculate eigenvectors and values. Recommended for very large matrices.
#' Set to NULL if you don't want to use this method
#' @param discrete_time Set to TRUE if movement model uses discrete time.
#' @param scale_by_density Set to TRUE to scale by the inverse of patch-specific long-term occupancy.
#' @param keep_imaginary Set to TRUE to retain complex values in calculation.
#' @param progress_bar Set to TRUE to view calculation progress.
#'
#' @returns TBD
#' @export
#'
#' @examples TBD
calculate_kinetic_distances = function(movement_matrix,
                                       d,
                                       T,
                                       sigma_val=1e-16,
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
    stop("If using a discrete-time movement matrix, the time scale T must be an integer")
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

#' Title
#'
#' @param cluster_type
#' @param landscape
#' @param out
#' @param min_dens
#' @param n_clust
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
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


calc_density <- function(disperse_mat, sigma=1e-16){
  dens = RSpectra::eigs(disperse_mat,k = 1,which = "LM",sigma = 1e-16)
  dens = Re(dens$vectors[,1])
  dens <- dens/sum(dens)
  dens
}
