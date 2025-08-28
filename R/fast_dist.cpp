#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::cx_colvec fastdist(arma::cx_mat Psi, 
                         arma::cx_mat inner) {
  // getting dimensions of the matrix
  int n = Psi.n_rows, d = Psi.n_cols;
  //there are n2 = n*(n-1)/2 possible unique pairs of n locations
  int n2 = n*(n-1)/2;
  int counter;
  //matrix of differences between values of the d left eigenvectors, in compact
  //form: each column i represents differences between the values of the ith
  //left eigenvector, excluding duplicated differences (so the first row is the 
  //difference between the eigenvector values for location 1 and 2, etc. 
  arma::cx_mat Psi_diff(n2, d);
  
  //creating the matrix of pairwise differences between mat values for each
  //unique pair of i and j for each dimension d. 
  for(int d_val = 0; d_val<d; d_val++){
    counter = 0;
    for(int i =0; i<(n-1); i++){
      for(int j = i+1; j<n; j++){
        //Calculates differences of the eigenvectors, scaled by the eigenvalues
        //for each column
        Psi_diff(counter, d_val) = (Psi(i,d_val) - Psi(j,d_val));
        counter++;
      }
    }
  }

  //creating empty matrices to store results in.
  arma::cx_mat dists(n2, d);
  arma::cx_colvec colsums(n2);
  
  //finding the distances between each pair of matrices
  for(int i = 0; i<d; i++){
    // this calculates the product of the matrix of differences with the ith
    //column of the Phi_inner matrix (the matrix of inner products of the 
    //right eigenvectors)
    colsums = Psi_diff*inner.col(i);
    
    //element-wise multiply the summed vector by the ith vector of differences
    dists.col(i) = colsums%Psi_diff.col(i);
   }
  
  //calculate the row-sum of the d contributions to the difference. 
  arma::cx_colvec dists_sums = sum(dists,1);
  return(dists_sums);
}