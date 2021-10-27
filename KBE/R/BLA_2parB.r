#' Bayes Linear Adjustment by 2 Parallel Known Boundaries
#'
#' @description Perform a Bayes linear adjustment utilising knowledge of function behaviour along
#' two parallel known boundaries in the input space.
#'
#' @param x points at which we want to update
#' @param K_d the dimensions which, when fixed at certain values, result in known boundary K.
#' @param L_d the dimensions which, when fixed at certain values, result in known boundary L.
#' @param xK the projection of x onto known boundary K
#' @param xL the projection of x onto known boundary L
#' @param xLK the projection of x first onto known boundary L and then known boundary K.
#' @param xK_d values the dimensions K must take for the function to be known
#' @param xL_d values the dimensions L must take for the function to be known
#' @param fxK function evaluated at x projected onto the boundary K.
#' @param fxL function evaluated at x projected onto the boundary L.
#' @param fxLK function evaluated at xLK.
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
#' L_d = c(2,3)
#' xK_d = 0
#' xL_d = c(1,1)
#' # If we are in a parallel setting, then xLK (projection of x first onto L and then K)
#' # is given as follows:
#' xLK_d <- xL_d
#' xLK_d[1:length(xK_d)] <- xK_d
#' # And LK_d (fixed values of coordinates for projections first onto L and then K)
#' # is just given by L_d.
#' LK_d <- L_d
#' fxK <- f_boundary( x = x, K_d = K_d, xK_d = xK_d, f = f )
#' fxL <- f_boundary( x = x, K_d = L_d, xK_d = xL_d, f = f )
#' fxLK <- f_boundary( x = x, K_d = LK_d, xK_d = xLK_d, f = f )
#' theta <- c( pi, pi/8, pi )
#' s2 <- 2
#' BA <- BLA_2parB( x = x, K_d = K_d, L_d = L_d, xK_d = xK_d, xL_d = xL_d,
#'                   fxK = fxK, fxL = fxL, fxLK = fxLK, theta = theta, s2 = s2 )
BLA_2parB <- function( x,
                       K_d, L_d = 0,
                       xK = NA, xL = NA, xLK = NA,
                       xK_d = NA, xL_d = NA,
                       fxK, fxL, fxLK,
                       E_fx = 0, E_fxK = 0, E_fxL = 0, E_fxLK = 0,
                       theta, s2 ){

  if( identical( K_d, L_d ) ){ LmK_d <- 0 }else{ LmK_d <- L_d[-(1:length(K_d))] }

  if( identical( xK, NA ) ){

    xK <- xL <- xLK <- x

    xLK[,L_d] <- xL[,L_d] <- rep( xL_d, times = rep( nrow(x), length(xL_d) ) )

    xLK[,K_d] <- xK[,K_d] <- rep( xK_d, times = rep( nrow(x), length(xK_d) ) )

  }

  # Correlation function between x's in dimension K_d and L_d-K_d.
  KXX_1k <- GaussianCF( X = x[,K_d], theta = theta[K_d] )
  KXX_kl <- GaussianCF( X = x[,LmK_d], theta = theta[LmK_d] )

  # Correlation function between x's in other dimensions.
  KXX_lp <- GaussianCF( X = x[,-L_d], theta = theta[-L_d] )

  # Correlation function for x in dimension kld between x and its boundary projections, xk and xl. r_1(a), r_1(b) etc.
  r_1k_a <- GaussianCF( X = x[,K_d], Y = xK[1,K_d,drop=FALSE], theta = theta[K_d] )
  r_1k_b <- GaussianCF( X = x[,K_d], Y = xL[1,K_d,drop=FALSE], theta = theta[K_d] )
  r_kl_b <- GaussianCF( X = x[,LmK_d], Y = xL[1,LmK_d,drop=FALSE], theta = theta[LmK_d] ) # this r_{k+1:l}(b) in the paper
  r_1k_LK <- GaussianCF( xK[1,K_d,drop=FALSE], xL[1,K_d,drop=FALSE], theta = theta[K_d] )[1,1]

  # Multiplication of each pair of correlation functions between x and projection, in known dimension. r_1(a)*r_1(a') etc.
  rr_1k_aa <- kronecker(r_1k_a, t(r_1k_a))
  rr_1k_bb <- kronecker(r_1k_b, t(r_1k_b))
  rr_kl_bb <- kronecker(r_kl_b, t(r_kl_b))

  # R functions
  R_1k_aa <- KXX_1k - rr_1k_aa
  R_1k_aLK <- r_1k_b - r_1k_a * r_1k_LK
  R_1k_LKLK <- 1 - r_1k_LK^2

  # Multiplication of each pair of R functions between x and projection, eg. R_1(a,LK)*R_1(LK,a') etc.
  RR_1k_aa <- kronecker(R_1k_aLK, t(R_1k_aLK))

  # R2 function
  R2_kl <- R_1k_aa * KXX_kl - ( RR_1k_aa * ( 1 / R_1k_LKLK ) * rr_kl_bb )

  # Updated expectation.
  ELK_fx <- E_fx +
    r_1k_a * ( fxK - E_fxK ) +
    R_1k_aLK * ( 1 / R_1k_LKLK ) * r_kl_b * ( ( fxL - E_fxL ) - r_1k_LK * ( fxLK - E_fxLK ) )

  # Updated Covariance
  CovLK_fx <- s2 * KXX_lp * R2_kl

  return( list( "EB_fx" = ELK_fx, "VarB_fx" = diag( CovLK_fx ), "CovB_fx" = CovLK_fx ) )

}
