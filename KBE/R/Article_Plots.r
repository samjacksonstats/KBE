#' Quick Generation of 3D Example Plots.
#'
#' @param f a toy function for which plots will be created.
#' @param ranges the ranges for the input parameters of the toy functions, given as a matrix.
#' @param K_d variable indices to be fixed for boundary K.
#' @param L_d variable indices to be fixed for boundary L.
#' @param M_d variable indices to be fixed for boundary M.
#' @param xK_d values at which those variables are fixed to.
#' @param xL_d values at which those variables are fixed to.
#' @param xM_d values at which those variables are fixed to.
#' @param fixed_dimension index of the variable that will be kept fixed for the plots, given as a vector of length 3.
#' @param fixed_value fixed value in the remaining dimension, given as a vector of length 3.
#' @param theta correlation length parameters.
#' @param s2 scalar variance parameter.
#' @param zlim_f plotting range for the z-values (model and mean prediction)
#' @param zlim_var plotting range for the variance predictions.
#' @param main_cube title for the cube plot, given as a list of length 3.
#'
#' @return nothing is returned.  Plots are produced.
#' @export
#'
#' @examples
#' # Specify the toy function - requires a 3D input vector.
#' f <- function( x ){
#'
#'  sin( x[1] / ( exp( x[2] ) ) ) + cos( x[3] )
#'
#' }
#'
#' # Specify ranges for the function parameters x1, x2, and x3.
#' ranges <- matrix(c( -2*pi, 2*pi,
#'                     -pi/4, pi/4,
#'                     -2*pi, 2*pi), ncol = 2, byrow = TRUE)
#'
#' # Specify the correlation length parameters and variance parameter for the example.
#' theta <- c( pi, pi/8, pi )
#' s2 <- 2
#'
#' # Specify the ranges for the function/mean and variance plots.
#' zlim_f <- c(-2.5, 2.5)
#' zlim_var <- c(0, 2)
#'
#' # Specify the boundary.
#' K_d <- c( 2, 3 )
#' xK_d <- c( 0, 0 )
#'
#' # Fixed dimensions.
#' fixed_dimension <- c( 2, 2, 1 )
#' fixed_value <- c( 0, -pi/8, -pi )
#'
#' # Set labels.
#' quotes <- list( bquote( x[2] == 0 ), bquote( x[2] == -pi/8 ), bquote( x[1] == -pi ) )
#'
#' # Run the function.
#' Article_Plots( f = f,
#'                ranges = ranges,
#'                K_d = K_d,
#'                xK_d = xK_d,
#'                fixed_dimension = fixed_dimension,
#'                fixed_value = fixed_value,
#'                theta = theta,
#'                s2 = s2,
#'                zlim_f = zlim_f,
#'                zlim_var = zlim_var,
#'                main_cube = quotes )
Article_Plots <- function( f,
                           ranges,
                           K_d,
                           L_d = NA,
                           M_d = NA,
                           xK_d,
                           xL_d = NA,
                           xM_d = NA,
                           fixed_dimension,
                           fixed_value,
                           theta = theta,
                           s2 = s2,
                           zlim_f = "assessed",
                           zlim_var = "assessed",
                           main_cube = list( "", "", "" ) ){

  # Constant vectors.
  legend_vec <- c( FALSE, FALSE, TRUE )

  # Set up the plot domain.
  KBE::plot_setup()

  # Loop over the three function plots.
  for( i in 1:3 ){

    KBE::PlotGen3dEx( f = f,
                      ranges = ranges,
                      K_d = K_d,
                      L_d = L_d,
                      M_d = M_d,
                      xK_d = xK_d,
                      xL_d = xL_d,
                      xM_d = xM_d,
                      fixed_dimension = fixed_dimension[i],
                      fixed_value = fixed_value[i],
                      theta = theta,
                      s2 = s2,
                      zlim_f = zlim_f,
                      zlim_var = zlim_var,
                      main_cube = main_cube[[i]],
                      legend = legend_vec[i] )

  }

}
