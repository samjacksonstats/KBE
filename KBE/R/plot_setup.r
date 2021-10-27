#' Plot Setup for Article.
#'
#' @description A simple function for setting up the plot domain to generate figures similar to those
#' presented for the 3D example shown in the article.
#'
#' @return Sets up plot domain.
#' @export
#'
#' @examples
#' plot_setup()
plot_setup <- function(){

  graphics::par( mfrow = c(5,5) )
  layout_matrix <- matrix( c(5:19, 0, 20:23, 0, 1:4), ncol = 5 )
  graphics::layout( layout_matrix, widths = c(5,5,5,2.4,1.4) )

  # Names for the plots, preferably given as expressions.
  plot_names <- list( expression( f(x) ), expression( mu( x ) ), expression( nu( x ) ), expression( s(x) ) )

  graphics::par( mar = c(0, 0.2, 0, 0.5) )
  for( j in 1:4 ){
    graphics::plot( 1, type = "n", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", bty = "n", xlim = c(0.5,5.5), ylim = c(0.5,5.5) )
    graphics::text( 3.5, 3, plot_names[[j]], srt = 0, cex = 1.8 )
  }

  graphics::par( mar = rep(1,4) )

}
