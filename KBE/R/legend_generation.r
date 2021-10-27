#' Legend Generation
#'
#' @description Generate the legend for the 3D example as shown in the article.
#'
#' @param colours colours for each level of z.
#' @param levels levels at which z should be divided into for the contour plot.
#'
#' @return nothing is returned.  Legend is plotted.
#' @export
#'
#' @examples
#' nlevels = 20
#' levels <- pretty( c(0,7), nlevels )
#' colours <- grDevices::colorRampPalette(
#'                c( "white", "cyan", "blue", "purple", "pink" ), space = "Lab" )
#' legend_generation( colours = colours, levels = levels )
legend_generation <- function( colours, levels ){

  graphics::plot( NA, ylim = range( levels ),
        xlim = c(0,1), xlab = "", ylab = "",
        frame = TRUE, axes = F, xaxs = "i", yaxs = "i" )
  graphics::axis( side = 4, cex.axis = 2, las = 2 )

  for( j in 1:( length( levels ) - 1 ) ){
    graphics::rect( 0, levels[j], 1, levels[j+1], col = colours( length( levels ) - 1 )[j] )
  }

}
