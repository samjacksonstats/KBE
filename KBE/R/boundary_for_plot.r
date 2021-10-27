#' Boundary Coordinate Generation
#'
#' @description A function for generating the coordinate matrices
#' (with 2 columns, and 4 rows for 2D boundaries and 2 rows for 1D boundaries)
#' required for the cube.
#'
#' @param fixed_dimension the dimensions which are fixed for the boundary
#  (this being a scalar for a 2D boundary or vector of length 2 for 1D boundary).
#' @param fixed_value the values for the coordinates which are fixed (of same length as fixed_dimension).
#' @param ranges ranges of the three variables, given as a 3 x 2 matrix.
#'
#' @return A matrix with 2 columns, and 4 rows for 2D boundaries and 2 rows for 1D boundaries
#' This gives the coordinates of the boundaries for plotting on the cube as shown in the top
#' row of the 3D example figures in the article.
#' @export
#'
#' @examples
#' ranges <- matrix(c( -2*pi, 2*pi,
#'                     -pi/4, pi/4,
#'                     -2*pi, 2*pi), ncol = 2, byrow = TRUE)
#' boundary_for_plot( fixed_dimension = 2,
#'                    fixed_value = 0,
#'                    ranges = ranges )
boundary_for_plot <- function( fixed_dimension,
                               fixed_value,
                               ranges ){

  # If there is only one fixed dimension, this means that we are plotting a two-dimensional boundary.
  if( length( fixed_dimension ) == 1 ){

    # Is the fixed dimension x1, x2 or x3?
    if( fixed_dimension == 1 ){
      sfv <- KBE::Scale( fixed_value, a = ranges[1,1], b = ranges[1,2], l = 1, u = 4 )[1,1] # scaled fixed value
      boundary <- matrix( c( sfv, 1,
                             sfv, 4,
                             sfv + 1, 5,
                             sfv + 1, 2 ), ncol = 2, byrow = TRUE )
    }

    if( fixed_dimension == 2 ){
      sfv <- KBE::Scale( fixed_value, a = ranges[2,1], b = ranges[2,2], l = 1, u = 2 )[1,1] # scaled fixed value
      boundary <- matrix( c( sfv, sfv,
                             sfv, sfv + 3,
                             sfv + 3, sfv + 3,
                             sfv + 3, sfv ), ncol = 2, byrow = TRUE )
    }

    if( fixed_dimension == 3 ){
      sfv <- KBE::Scale( fixed_value, a = ranges[3,1], b = ranges[3,2], l = 1, u = 4 )[1,1] # scaled fixed value
      boundary <- matrix(c( 1, sfv,
                            4, sfv,
                            5, sfv + 1,
                            2, sfv + 1 ), ncol = 2, byrow = TRUE )
    }

  }

  # If there are two fixed dimensions, this means that we are plotting a one-dimensional boundary.
  if( length( fixed_dimension ) == 2 ){

    # Are the fixed dimensions (1,2), (1,3) or (2,3)?
    if( identical( fixed_dimension, c(1,2) ) ){
      sfv <- KBE::Scale( matrix( fixed_value, nrow = 1 ), a = ranges[1:2,1], b = ranges[1:2,2], l = c(1,0), u = c(4,1) )
      boundary <- matrix( c( sfv[1] + sfv[2], 1 + sfv[2],
                             sfv[1] + sfv[2], 4 + sfv[2] ), ncol = 2, byrow=TRUE )
    }

    if( identical( fixed_dimension, c(1,3) ) ){
      sfv <- KBE::Scale( matrix( fixed_value, nrow = 1 ), a=ranges[c(1,3),1], b = ranges[c(1,3),2], l=c(1,1), u=c(4,4) )
      boundary <- matrix(c( sfv[1], sfv[2],
                            sfv[1] + 1, sfv[2] + 1), ncol = 2, byrow=TRUE)
    }

    if( identical( fixed_dimension, c(2,3) ) ){
      sfv <- KBE::Scale( matrix( fixed_value, nrow = 1 ), a = ranges[2:3,1], b = ranges[2:3,2], l = c(0,1), u = c(1,4) )
      boundary <- matrix( c( sfv[1] + 1, sfv[1] + sfv[2],
                             sfv[1] + 4, sfv[1] + sfv[2] ), ncol = 2, byrow=TRUE )
    }

  }

  return( boundary )

}
