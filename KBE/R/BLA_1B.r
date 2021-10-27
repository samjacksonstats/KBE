#' Bayes Linear Adjustment by a Single Known Boundary
#'
#' @description Perform a Bayes linear adjustment utilising knowledge of function behaviour along
#' a single boundary in the input space.
#'
#' @param x points at which we want to update
#' @param K_d the dimensions which, when fixed at certain values, result in known boundaries.
#' @param xK the projection of x onto known boundary K
#' @param xK_d values the dimensions K must take for the function to be known
#' @param fxK function evaluated at x projected onto the boundary K.
#' @param E_fx prior expectation for the function f(x)
#' @param E_fxK prior expectation for f(x^K)
#' @param theta vector of correlation length parameter values
#' @param s2 scalar variance parameter value.
#'
#' @return
#' \item{EB_fx}{Expected value of f(x) adjusted by knowledge of function behaviour along K.}
#' \item{VarB_fx}{Variance of f(x) adjusted by knowledge of function behaviour along K.}
#' \item{CovB_fx}{Covariance of f(x) adjusted by knowledge of function behaviour along K.}
#' @export
#'
#' @examples
#' # Toy function
#' f <- function( x ){
#'
#' sin( x[1] / ( exp( x[2] ) ) ) + cos( x[3] )
#'
#' }
#' x <- matrix( runif( 12 ), ncol = 3 )
#' K_d = 2
#' xK_d = 0
#' fxK <- f_boundary( x = x, K_d = K_d, xK_d = xK_d, f = f )
#' theta <- c( pi, pi/8, pi )
#' s2 <- 2
#' BA <- BLA_1B( x = x, K_d = K_d, xK_d = xK_d, fxK = fxK, theta = theta, s2 = s2 )
BLA_1B <- function( x, K_d, xK = NA, xK_d = NA, fxK, E_fx = 0, E_fxK = 0, theta, s2 ){

  # Set the matrix of projections of x onto the known boudnary, if not specified.
  if( identical( xK, NA ) ){

    xK <- x

    xK[,K_d] <- rep( xK_d, times = rep( nrow(x), length(xK_d) ) )

  }

  # Correlation function between x's in dimension K_d which have to take fixed values.
  KXX_1k <- KBE::GaussianCF( X = x[,K_d], theta = theta[K_d] )

  # Correlation function between x's in other dimensions.
  KXX_kp <- KBE::GaussianCF( X = x[,-K_d], theta = theta[-K_d] )

  # Correlation function for x in dimension K_d between x and its boundary projection, xk.
  r_1k_a <- KBE::GaussianCF( X = x[,K_d], Y = xK[1,K_d,drop=FALSE], theta = theta[K_d] )

  # Multiplication of each pair of correlation functions between x and projection, in known dimension.
  rr_1k_aa <- r_1k_a %*% t(r_1k_a)

  # Difference between KXX_1k and rr_1k_aa: R_1(a,a') = r_1(a-a') - r_1(a) * r_1(a').
  R_1k_aa <- KXX_1k - rr_1k_aa

  # Adjusted expectation.
  EK_fx <- E_fx + r_1k_a * ( fxK - E_fxK )

  # Adjusted Covariance.
  CovK_fx <- s2 * R_1k_aa * KXX_kp

  # Return list of required objects.
  return( list( "EB_fx" = EK_fx, "VarB_fx" = diag( CovK_fx ), "CovB_fx" = CovK_fx ) )

}
