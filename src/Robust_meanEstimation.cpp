// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
using namespace arma;

//' @title High-Dimensional Robust Mean Estimation via Gradient Descent using Rcpp
//' @description High-Dimensional Robust Mean Estimation via Gradient Descent using Rcpp
//' @param matrix Design matrix over \eqn{R^{p \times n}}, \code{n} is the number of random vectors, \code{p} is the dimension of random vectors
//' @param epsilon Estimated proportion of outliers in a dataset
//' @param max_iter The maximum number of iterations of the algorithm
//' @param tol The convergence coefficient of the algorithm
//' @return A weight vector for a weighted mean estimate
//' @examples
//' \dontrun{
//' random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
//' w <- projected_gradient_descent(random_matrix_normal, epsilon)
//' mu_w <- random_matrix_normal %*% w
//' print("The weighted mean estimate:")
//' }
//' @export
// [[Rcpp::export]]
arma::vec projected_gradient_descent(arma::mat matrix, double epsilon, int max_iter = 1000, double tol = 1e-5) {
  int d = matrix.n_rows;  // the dimension of a dataset
  int N = matrix.n_cols;  // the number of sample
  double eta = 1.0 / (2 * std::pow(arma::norm(matrix, "inf"), 2));  // step size
  arma::vec w(N, fill::ones);  // The initialization weight vector is evenly distributed
  w /= N;
  double delta = 2 * epsilon * N;  // The maximum allowable weight for each sample
  for (int iter = 0; iter < max_iter; ++iter) {
    // compute weighted mean
    arma::vec mu_w = matrix * w;
    // compute weighted covariance matrix
    arma::mat matrix_centered = matrix.each_col() - mu_w;
    arma::mat Sigma_w = matrix_centered * diagmat(w) * matrix_centered.t() - mu_w * mu_w.t();
    // compute gradient
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Sigma_w);
    arma::vec u = eigvec.col(d - 1);  // Maximum eigenvector
    arma::vec grad(N);
    for (int i = 0; i < N; ++i) {
      grad[i] = arma::dot(matrix_centered.col(i), u) * arma::dot(matrix_centered.col(i), u) - 2 * arma::dot(matrix_centered.col(i), u) * arma::dot(mu_w, u);
    }
    // Steps for gradient descent
    arma::vec w_new = w - eta * grad;
    // Projection to a feasible set
    w_new = arma::clamp(w_new, 0.0, delta);
    w_new /= arma::sum(w_new);
    // Check for convergence
    if (arma::norm(w_new - w, 2) < tol) {
      break;
    }
    w = w_new;
  }
  return w;
}