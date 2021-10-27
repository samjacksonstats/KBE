#' Bayes Linear Adjustment by 2 Perpendicular Known Boundaries
#'
#' @description Perform a Bayes linear adjustment utilising knowledge of function behaviour along
#' two perpendicular known boundaries in the input space.
#'
#' @param x points at which we want to update
#' @param K_d the dimensions which, when fixed at certain values, result in known boundary K.
#' @param L_d the dimensions which, when fixed at certain values, result in known boundary L.
#' @param xK the projection of x onto known boundary K
#' @param xL the projection of x onto known boundary L
#' @param xLK the projection of x onto the intersection of known boundaries K and L.
#' @param xK_d values the dimensions K must take for the function to be known
#' @param xL_d values the dimensions L must take for the function to be known
#' @param fxK function evaluated at x projected onto the boundary K.
#' @param fxL function evaluated at x projected onto the boundary L.
#' @param fxLK function evaluated at x projected onto the intersection of boundaries K and L.
#' @param E_fx prior expectation for the function f(x)
#' @param E_fxK prior expectation for f(x^K)
#' @param E_fxL prior expectation for f(x^L)
#' @param E_fxLK prior expectation for f(x^LK)
#' @param theta vector of correlation length parameter values.
#' @param s2 scalar variance parameter value.
#'
#' @return
#' \item{EB_fx}{Expected value of f(x) adjusted by knowledge of function behaviour along K and L.}
#' \item{VarB_fx}{Variance of f(x) adjusted by knowledge of function behaviour along K and L.}
#' \item{CovB_fx}{Covariance of f(x) adjusted by knowledge of function behaviour along K and L.}
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
#' L_d = 1
#' xK_d = 0
#' xL_d = 0
#' fxK <- f_boundary( x = x, K_d = K_d, xK_d = xK_d, f = f )
#' fxL <- f_boundary( x = x, K_d = L_d, xK_d = xL_d, f = f )
#' fxLK <- f_boundary( x = x, K_d = c(K_d, L_d), xK_d = c(xK_d, xL_d), f = f )
#' theta <- c( pi, pi/8, pi )
#' s2 <- 2
#' BA <- BLA_2perpB( x = x, K_d = K_d, L_d = L_d, xK_d = xK_d, xL_d = xL_d,
#'                   fxK = fxK, fxL = fxL, fxLK = fxLK, theta = theta, s2 = s2 )
BLA_2perpB <- function(x, K_d, L_d, xK = NA, xL = NA, xLK = NA, xK_d = NA, xL_d = NA, fxK, fxL, fxLK, E_fx=0, E_fxK=0, E_fxL=0, E_fxLK=0, theta, s2){

  if( identical( xK, NA ) ){

    xK <- xL <- xLK <- x

    xLK[,K_d] <- xK[,K_d] <- rep( xK_d, times = rep( nrow(x), length(xK_d) ) )

    xLK[,L_d] <- xL[,L_d] <- rep( xL_d, times = rep( nrow(x), length(xL_d) ) )

  }

  # Correlation function between x's in dimension kd and ld.
  KXX_1k <- KBE::GaussianCF( X = x[,K_d], theta = theta[K_d] )
  KXX_kl <- KBE::GaussianCF( X = x[,L_d], theta = theta[L_d] )

  # Correlation function between x's in other dimensions.
  if( ncol( x ) == sum( length( K_d ) + length( L_d ) ) ){
    KXX_lp <- matrix( 1, nrow = nrow( x ), ncol = nrow( x ) )
  }else{
    KXX_lp <- KBE::GaussianCF( X = x[,-c(K_d,L_d)], theta = theta[-c(K_d,L_d)] )
  }

  # Correlation function for x in dimension kd between x and its boundary projection, xk. r_1(a), r_2(b) etc.
  r_1k_a <- KBE::GaussianCF( X = x[,K_d], Y = xK[1,K_d,drop=FALSE], theta = theta[K_d] )
  r_kl_b <- KBE::GaussianCF( X = x[,L_d], Y = xL[1,L_d,drop=FALSE], theta = theta[L_d] )

  # Multiplication of each pair of correlation functions between x and projection, in known dimension. r_1(a)*r_1(a') etc.
  rr_1k_aa <- r_1k_a %*% t( r_1k_a )
  rr_kl_bb <- r_kl_b %*% t( r_kl_b )

  # Difference between KXXkd and rramat.  R_1(a,a') = r_1(a-a') - r_1(a)*r_1(a') etc.
  R_1k_aa <- KXX_1k - rr_1k_aa
  R_kl_bb <- KXX_kl - rr_kl_bb

  # Updated expectation.
  ELK_fx <- E_fx + r_1k_a * ( fxK - E_fxK ) +
    r_kl_b * ( fxL - E_fxL ) +
    -r_1k_a * r_kl_b * ( fxLK - E_fxLK )

  # Updated Covariance
  CovLK_fx <- s2 * R_1k_aa * R_kl_bb * KXX_lp

  return( list( "EB_fx" = ELK_fx, "VarB_fx" = diag( CovLK_fx ), "CovB_fx" = CovLK_fx ) )

}
