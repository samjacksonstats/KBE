#' Draw a Cube
#'
#' @description Specific function for the diagrams of a cube illustrating which boundaries are known and
#' which cross-section of the input space is being emulated.
#'
#' @param lwd line width for the edges of the cube.
#' @param lty2 line type for the "hidden" edges of the cube.
#' @param cex scale size for the labels.
#' @param col_line_width line width for the coloured lines representing 1D boundaries.
#' @param coloured_hp a list of vectors, where each vector represents the hyperplanes to be plotted:
#                these being of length 8 (4 coordinate pairs) for 2D boundaries
#                and of length 4 (2 coordinate pairs) for 1D boundaries.
#' @param col_hp the colour of the hyperplanes given above.
#' @param density_col the density of the fill of the 2D hyperplanes (note that the length of this vector
#               still needs to be equal to the length of coloured_hp even though the elements
#               corresponding to 1D boundaries won't be used.
#' @param main title of the cube plot.
#' @param cex.main size of the title.
#' @param main.line number of lines outwards from the plot edge to plot the title.
#'
#' @return nothing is returned.  Cube is plotted.
#' @export
#'
#' @examples
#' ranges <- matrix(c( -2*pi, 2*pi,
#'                     -pi/4, pi/4,
#'                     -2*pi, 2*pi), ncol = 2, byrow = TRUE)
#' boundary_for_plot <- boundary_for_plot( fixed_dimension = 2,
#'                                         fixed_value = 0,
#'                                         ranges = ranges )
#' draw_cube( coloured_hp = list( boundary_for_plot ) )
draw_cube <- function( lwd = 2,
                       lty2 = 2,
                       cex = 1.8,
                       col_line_width = lwd,
                       coloured_hp = list(),
                       col_hp = c( "green", "red", "blue", "pink" ),
                       density_col = rep( 0.7, length( coloured_hp ) ),
                       main = "",
                       cex.main = 1,
                       main.line = 1 ){

  # Empty plotting device.
  graphics::plot( 1, type = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i",
                  xlab = "", ylab = "", bty = "n", xlim = c(0.3,5.7), ylim = c(0.3,5.7) )
  graphics::title( main = main, cex.main = cex.main, line = main.line )

  # Number of hyperplanes
  n_hyp <- length( coloured_hp )

  # Plot the 2D Boundaries, if required.
  for( k in 1:n_hyp ){
    if( length( coloured_hp[[k]] ) == 8 ){
      graphics::polygon( coloured_hp[[k]], col = grDevices::adjustcolor( col_hp[k], alpha = density_col[k] ), border = "NA" )
    }
  }

  # Coordinates of the cube outline as a 2D shape.
  coordinates <- matrix( c( 1, 1,
                            1, 4,
                            2, 5,
                            5, 5,
                            5, 2,
                            4, 1 ), ncol = 2, byrow = TRUE )

  # Plot shape outline as given by the coordinates above.
  graphics::polygon( coordinates, lwd = lwd )

  # Additional lines to fill out the cube.
  graphics::lines( c(1,4), c(4,4), lwd = lwd )
  graphics::lines( c(4,4), c(4,1), lwd = lwd )
  graphics::lines( c(4,5), c(4,5), lwd = lwd )
  graphics::lines( c(1,2), c(1,2), lty = lty2, lwd = lwd )
  graphics::lines( c(2,5), c(2,2), lty = lty2, lwd = lwd )
  graphics::lines( c(2,2), c(2,5), lty = lty2, lwd = lwd )

  # Add dimension text to the cube.
  graphics::text( c( 2.4, 1.4, 0.68 ),
                  c( 0.7, 1.95, 3 ),
                  c( expression( x[1] ), expression( x[2] ), expression( x[3] ) ),
                  cex = cex )

  # Plot 1D boundaries, if required.
  for( k in 1:n_hyp ){
    if( length( coloured_hp[[k]] ) == 4 ){ graphics::lines( coloured_hp[[k]], col = col_hp[k],  lwd = col_line_width ) }
  }

}
