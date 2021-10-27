#' Evaluate f along a boundary.
#'
#' @param x set of points, given as a matrix.
#' @param K_d Variable indices to be fixed.
#' @param xK_d Values at which those variables are fixed to.
#' @param f a toy function which is to be evaluated.
#'
#' @return Value of f(x) at the projections of x projected onto boundary K_d.
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
#' f_boundary( x = x,
#'             K_d = 2,
#'             xK_d = 0,
#'             f = f )
f_boundary <- function( x, K_d, xK_d, f ){

  x[,K_d] <- rep( xK_d, times = rep( nrow( x ), length( xK_d ) ) )

  apply( x, 1, f )

}
