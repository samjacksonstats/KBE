#' Generate Contour Plot
#'
#' @description User-friendly wrapper for generating the contour plots as of the 3D example shown in the text.
#'
#' @param x x coordinates.
#' @param y y coordinates.
#' @param z matrix of values to be plotted, with rows assumed to correspond to increasing values of x (from top to bottom)
#     and columns to y (from left to right).
#' @param levels levels at which z should be divided into for the contour plot.
#' @param colours colours for each level of z.
#'
#' @return nothing is returned.  Contour plot is generated.
#' @export
#'
#' @examples
#' x = seq( from = 0 , to = 1, by = 0.01 )
#' y = seq( from = 0 , to = 1, by = 0.01 )
#' eg <- expand.grid( x, y )
#' z <- matrix( apply( eg, 1, FUN = function(x)( cos(x[1]) + tan(x[2]) )^2 ), nrow = 101 )
#' nlevels = 20
#' levels <- pretty( c(0,7), nlevels )
#' colours <- grDevices::colorRampPalette(
#'                c( "white", "cyan", "blue", "purple", "pink" ), space = "Lab" )
#' contour_plot( x = x, y = y, z = z, levels = levels, colours = colours )
contour_plot <- function( x, y, z, levels, colours ){

  # Set up an empty plot of the right dimensions.
  graphics::plot( NA,
                  xlim = range( x ),
                  ylim = range( y ),
                  xlab = "",
                  ylab = "",
                  frame = TRUE,
                  axes = F,
                  xaxs = "i",
                  yaxs = "i",
                  main = "" )


  graphics::.filled.contour( x = x,
                             y = y,
                             z = z,
                             levels = levels,
                             col = colours( length( levels ) - 1 ) )

}
