#' Bayes Linear Adjustment by 3 Perpendicular Known Boundaries
#'
#' @description Perform a Bayes linear adjustment utilising knowledge of function behaviour along
#' three perpendicular known boundaries in the input space.
#'
#' @param x points at which we want to update
#' @param K_d the dimensions which, when fixed at certain values, result in known boundary K.
#' @param L_d the dimensions which, when fixed at certain values, result in known boundary L.
#' @param M_d the dimensions which, when fixed at certain values, result in known boundary M.
#' @param xK the projection of x onto known boundary K
#' @param xL the projection of x onto known boundary L
#' @param xM the projection of x onto known boundary M
#' @param xLK the projection of x onto the intersection of known boundaries K and L.
#' @param xMK the projection of x onto the intersection of known boundaries K and M.
#' @param xML the projection of x onto the intersection of known boundaries L and M.
#' @param xMLK the projection of x onto the intersection of known boundaries K, L and M.
#' @param xK_d values the dimensions K must take for the function to be known
#' @param xL_d values the dimensions L must take for the function to be known
#' @param xM_d values the dimensions M must take for the function to be known
#' @param fxK function evaluated at x projected onto the boundary K.
#' @param fxL function evaluated at x projected onto the boundary L.
#' @param fxM function evaluated at x projected onto the boundary M.
#' @param fxLK function evaluated at x projected onto the intersection of boundaries K and L.
#' @param fxMK function evaluated at x projected onto the intersection of boundaries K and M.
#' @param fxML function evaluated at x projected onto the intersection of boundaries L and M.
#' @param fxMLK function evaluated at x projected onto the intersection of boundaries K, L and M.
#' @param E_fx prior expectation for the function f(x)
#' @param E_fxK prior expectation for f(x^K)
#' @param E_fxL prior expectation for f(x^L)
#' @param E_fxM prior expectation for f(x^M)
#' @param E_fxLK prior expectation for f(x^LK)
#' @param E_fxMK prior expectation for f(x^MK)
#' @param E_fxML prior expectation for f(x^ML)
#' @param E_fxMLK prior expectation for f(x^MLK)
#' @param theta vector of correlation length parameter values.
#' @param s2 scalar variance parameter value.
#'
#' @return
#' \item{EB_fx}{Expected value of f(x) adjusted by knowledge of function behaviour along K, L and M.}
#' \item{VarB_fx}{Variance of f(x) adjusted by knowledge of function behaviour along K, L and M.}
#' \item{CovB_fx}{Covariance of f(x) adjusted by knowledge of function behaviour along K, L and M.}
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
#' M_d = 3
#' xK_d = 0
#' xL_d = 0
#' xM_d = 0
#' fxK <- f_boundary( x = x, K_d = K_d, xK_d = xK_d, f = f )
#' fxL <- f_boundary( x = x, K_d = L_d, xK_d = xL_d, f = f )
#' fxM <- f_boundary( x = x, K_d = M_d, xK_d = xM_d, f = f )
#' fxLK <- f_boundary( x = x, K_d = c(K_d, L_d), xK_d = c(xK_d, xL_d), f = f )
#' fxMK <- f_boundary( x = x, K_d = c(K_d, M_d), xK_d = c(xK_d, xM_d), f = f )
#' fxML <- f_boundary( x = x, K_d = c(L_d, M_d), xK_d = c(xL_d, xM_d), f = f )
#' fxMLK <- f_boundary( x = x, K_d = c(K_d, L_d, M_d), xK_d = c(xK_d, xL_d, xM_d), f = f )
#' theta <- c( pi, pi/8, pi )
#' s2 <- 2
#' BA <- BLA_3perpB( x = x, K_d = K_d, L_d = L_d, M_d = M_d,
#'                   xK_d = xK_d, xL_d = xL_d, xM_d = xM_d,
#'                   fxK = fxK, fxL = fxL, fxM = fxM,
#'                   fxLK = fxLK, fxMK = fxMK, fxML = fxML, fxMLK = fxMLK,
#'                   theta = theta, s2 = s2 )
BLA_3perpB <- function(x,
                       K_d, L_d, M_d,
                       xK = NA, xL = NA, xM = NA, xLK = NA, xMK = NA, xML = NA, xMLK = NA,
                       xK_d = NA, xL_d = NA, xM_d = NA,
                       fxK, fxL, fxM, fxLK, fxMK, fxML, fxMLK,
                       E_fx=0, E_fxK=0, E_fxL=0,E_fxM=0, E_fxLK=0, E_fxMK=0, E_fxML=0, E_fxMLK=0,
                       theta, s2){

  if( identical( xK, NA ) ){

    xK <- xL <- xM <- xLK <- xMK <- xML <- xMLK <- x

    xMLK[,K_d] <- xMK[,K_d] <- xLK[,K_d] <- xK[,K_d] <- rep( xK_d, times = rep( nrow(x), length(xK_d) ) )

    xMLK[,L_d] <- xML[,L_d] <- xLK[,L_d] <- xL[,L_d] <- rep( xL_d, times = rep( nrow(x), length(xL_d) ) )

    xMLK[,M_d] <- xML[,M_d] <- xMK[,M_d] <- xM[,M_d] <- rep( xM_d, times = rep( nrow(x), length(xM_d) ) )

  }

  # Correlation function between x's in dimension kd and ld.
  KXX_1k <- KBE::GaussianCF(X = x[,K_d], theta = theta[K_d])
  KXX_kl <- KBE::GaussianCF(X = x[,L_d], theta = theta[L_d])
  KXX_lm <- KBE::GaussianCF(X = x[,M_d], theta = theta[M_d])

  # Correlation function between x's in other dimensions.
  if(ncol(x) == sum( length( K_d ) + length( L_d ) + length( M_d ) ) ){
    KXX_mp <- matrix(1, nrow=nrow(x), ncol=nrow(x))
  }else{
    KXX_mp <- KBE::GaussianCF(X = x[,-c(K_d,L_d,M_d)], theta = theta[-c(K_d,L_d,M_d)])
  }

  # Correlation function for x in dimension kd between x and its boundary projection, xk. r_1(a), r_2(b) etc.
  r_1k_a <- KBE::GaussianCF( X = x[,K_d], Y = xK[1,K_d,drop=FALSE], theta = theta[K_d] )
  r_kl_b <- KBE::GaussianCF( X = x[,L_d], Y = xL[1,L_d,drop=FALSE], theta = theta[L_d] )
  r_lm_c <- KBE::GaussianCF( X = x[,M_d], Y = xM[1,M_d,drop=FALSE], theta = theta[M_d] )

  # Multiplication of each pair of correlation functions between x and projection, in known dimension. r_1(a)*r_1(a') etc.
  rr_1k_aa <- r_1k_a %*% t( r_1k_a )
  rr_kl_bb <- r_kl_b %*% t( r_kl_b )
  rr_lm_cc <- r_lm_c %*% t( r_lm_c )

  # Difference between KXXkd and rramat.  R_1(a,a') = r_1(a-a') - r_1(a)*r_1(a') etc.
  R_1k_aa <- KXX_1k - rr_1k_aa
  R_kl_bb <- KXX_kl - rr_kl_bb
  R_lm_cc <- KXX_lm - rr_lm_cc

  # Updated expectation.
  EMLK_fx <- E_fx + r_1k_a * ( fxK - E_fxK ) +
    r_kl_b * ( fxL - E_fxL ) +
    r_lm_c * ( fxM - E_fxM ) +
    -r_1k_a * r_kl_b * ( fxLK - E_fxLK ) +
    -r_1k_a * r_lm_c * ( fxMK - E_fxMK ) +
    -r_kl_b * r_lm_c * ( fxML - E_fxML ) +
    r_1k_a * r_kl_b * r_lm_c * ( fxMLK - E_fxMLK )

  # Updated Covariance
  CovMLK_fx <- s2 * R_1k_aa * R_kl_bb * R_lm_cc * KXX_mp

  return(list( "EB_fx" = EMLK_fx, "VarB_fx" = diag( CovMLK_fx ), "CovB_fx" = CovMLK_fx ) )

}
