#' 3D Example Plot Generation
#'
#' @description  Generic Generation of Plots for 3D Toy Example.
#'
#' @param f a toy function for which plots will be created.
#' @param ranges  the ranges for the input parameters of the toy functions, given as a matrix.
#' @param K_d variable indices to be fixed for boundary K.
#' @param L_d variable indices to be fixed for boundary L.
#' @param M_d variable indices to be fixed for boundary M.
#' @param xK_d values at which those variables are fixed to.
#' @param xL_d values at which those variables are fixed to.
#' @param xM_d values at which those variables are fixed to.
#' @param fixed_dimension index of the variable that will be kept fixed for the plots.
#' @param fixed_value fixed value in the remaining dimension.
#' @param theta correlation length parameters.
#' @param s2 scalar variance parameter.
#' @param grid_length number of grid points along each dimension with which to represent the plotted surface.
#' @param lwd line width for cube edges.
#' @param lty2 line type for background lines
#' @param cex relative size of plot
#' @param col_line_width width of coloured boundary lines
#' @param zlim_f plotting range for the z-values (model and mean prediction)
#' @param zlim_var plotting range for the variance predictions.
#' @param zlim_diag plotting range for the diagnostic plot.
#' @param legend Add a legend to the side of the plot?
#' @param main_cube title for the cube plot.
#' @param cex.main_cube font size of the cube plot title.
#' @param main.line line value for the title function.
#'
#' @return nothing is returned.  Plots are generated.
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
#' plot_setup()
#' PlotGen3dEx( f = f,
#'              ranges = ranges,
#'              K_d = c(2,3),
#'              L_d = c(2,3),
#'              M_d = 1,
#'              xK_d = c(0,0),
#'              xL_d = c(0,-pi),
#'              xM_d = 0,
#'              fixed_dimension = 2,
#'              fixed_value = 0,
#'              theta = theta,
#'              s2 = s2,
#'              zlim_f = zlim_f,
#'              zlim_var = zlim_var,
#'              main_cube = bquote( x[2] == 0 ),
#'              legend = FALSE )
PlotGen3dEx <- function( f,
                         ranges,
                         K_d,
                         L_d = NA,
                         M_d = NA,
                         xK_d,
                         xL_d = NA,
                         xM_d = NA,
                         fixed_dimension,
                         fixed_value,
                         theta = rep( pi/2, 3 ),
                         s2 = 1,
                         grid_length = 50,
                         lwd = 2,
                         lty2 = 2,
                         cex = 1.8,
                         col_line_width = 3,
                         zlim_f = "assessed",
                         zlim_var = "assessed",
                         zlim_diag = c(-4.25, 4.25),
                         legend = FALSE,
                         main_cube = "",
                         cex.main_cube = 1.8,
                         main.line = 0.2
){

  ## Diagnostic Slice Specification ##

  # Ranges over which the non-fixed variables of the diagnostic slice vary over.
  ranges_NF <- ranges[-fixed_dimension, ]

  # The 2D diagnostic slice is therefore defined by the coordinates of diagnostic_slice below...
  x1r <- seq( from = ranges_NF[1,1], to = ranges_NF[1,2], length = grid_length )
  x2r <- seq( from = ranges_NF[2,1], to = ranges_NF[2,2], length = grid_length )
  xp <- expand.grid( x1r, x2r )
  diagnostic_slice <- matrix( fixed_value, nrow = nrow( xp ), ncol = 3 )
  diagnostic_slice[,-fixed_dimension] <- as.matrix( xp )

  ## Evaluate f across diagnostic slice ##

  # Evaluate f at each point of diagnostic_slice for diagnostic purposes.
  fxu <- apply( diagnostic_slice, 1, f )

  # Put the model output into a matrix such that the rows correspond to increasing values of the first non-fixed variable,
  # and the columns correspond to increasing values of the second non-fixed variable.
  mdxv <- matrix( fxu, nrow = grid_length, ncol = grid_length, byrow = FALSE )

  ## Specification of Multiple Boundary Projections ##

  # If a second boundary L has been specified...
  if( identical( L_d, NA ) == FALSE ){

    # Are K and L perpendicular or parallel to each other?
    if( identical( K_d, L_d[1:length( K_d )] ) ){ LK <- "parallel" }else{ LK <- "perpendicular" }

    # If LK are parallel to each other...
    if( LK == "parallel" ){

      # If we are in a parallel setting, then xLK (projection of x first onto L and then K) is given as follows.
      xLK_d <- xL_d
      xLK_d[1:length(xK_d)] <- xK_d

      # And LK_d (fixed values of coordinates for projections first onto L and then K) is just given by L_d.
      LK_d <- L_d

      # Otherwise, if LK are perpendicular to each other...
    }else{

      LK_d <- c( K_d, L_d )
      xLK_d <- c( xK_d, xL_d )

    }

  }

  # If a third boundary M has been specified...
  if( identical( M_d, NA ) == FALSE ){

    # M is perpendicular either way...
    MK_d <- c( K_d, M_d )
    ML_d <- c( L_d, M_d )
    MLK_d <- c( LK_d, M_d )
    xMK_d <- c( xK_d, xM_d )
    xML_d <- c( xL_d, xM_d )
    xMLK_d <- c( xLK_d, xM_d )

  }

  ## Evaluate f at Boundaries ##

  # Use f_boundary to evaluate f at the projections of the grid of points diagnostic_slice above
  # (representing the diagnostic slice) onto the boundary defined by fixing variables K_d at values xK_d.
  fxuK <- KBE::f_boundary( x = diagnostic_slice, K_d = K_d, xK_d = xK_d, f = f )

  if( identical( L_d, NA ) == FALSE ){
    fxuL <- KBE::f_boundary( x = diagnostic_slice, K_d = L_d, xK_d = xL_d, f = f )
    fxuLK <- KBE::f_boundary( x = diagnostic_slice, K_d = LK_d,  xK_d = xLK_d, f = f )
  }

  if( identical( M_d, NA ) == FALSE ){
    fxuM <- KBE::f_boundary( x = diagnostic_slice, K_d = M_d, xK_d = xM_d, f = f )
    fxuMK <- KBE::f_boundary( x = diagnostic_slice, K_d = MK_d, xK_d = xMK_d, f = f )
    fxuML <- KBE::f_boundary( x = diagnostic_slice, K_d = ML_d, xK_d = xML_d, f = f )
    fxuMLK <- KBE::f_boundary( x = diagnostic_slice, K_d = MLK_d, xK_d = xMLK_d, f = f )
  }

  ## BL Boundary Adjustment ##

  # Depending on how many boundaries we have, and which combination they are (perpendicular or parallel),
  # we use a different KBE calculation function.

  if( identical( M_d, NA ) == FALSE ){

    if( identical( LK, "perpendicular" ) ){

      BL_boundary_adj <- KBE::BLA_3perpB( x = diagnostic_slice,
                                          K_d = K_d, L_d = L_d, M_d = M_d,
                                          xK_d = xK_d, xL_d = xL_d, xM_d = xM_d,
                                          fxK = fxuK, fxL = fxuL, fxM = fxuM,
                                          fxLK = fxuLK, fxMK = fxuMK, fxML = fxuML,
                                          fxMLK = fxuMLK,
                                          theta = theta, s2 = s2 )

    }else{

      BL_boundary_adj <- KBE::BLA_3B( x = diagnostic_slice,
                                      K_d = K_d, L_d = L_d, M_d = M_d,
                                      xK_d = xK_d, xL_d = xL_d, xM_d = xM_d,
                                      fxK = fxuK, fxL = fxuL, fxM = fxuM,
                                      fxLK = fxuLK, fxMK = fxuMK, fxML = fxuML,
                                      fxMLK = fxuMLK,
                                      theta = theta, s2 = s2 )

    }

  }else{

    if( identical( L_d, NA ) == FALSE ){

      if( identical( LK, "perpendicular" ) ){

        BL_boundary_adj <- KBE::BLA_2perpB( x = diagnostic_slice,
                                            K_d = K_d, L_d = L_d,
                                            xK_d = xK_d, xL_d = xL_d,
                                            fxK = fxuK, fxL = fxuL, fxLK = fxuLK,
                                            theta = theta, s2 = s2 )

      }else{

        BL_boundary_adj <- KBE::BLA_2parB( x = diagnostic_slice,
                                           K_d = K_d, L_d = L_d,
                                           xK_d = xK_d, xL_d = xL_d,
                                           fxK = fxuK, fxL = fxuL, fxLK = fxuLK,
                                           theta = theta, s2 = s2 )

      }

    }else{

      BL_boundary_adj <- KBE::BLA_1B( x = diagnostic_slice,
                                      K_d = K_d,
                                      xK_d = xK_d,
                                      fxK = fxuK,
                                      theta = theta, s2 = s2 )

    }

  }

  # Put the results of the adjustment above into corresponding matrices; one for expectation and one for variance.
  EB_fx <- matrix( BL_boundary_adj$EB_fx, nrow = grid_length, ncol = grid_length, byrow = FALSE )
  VarB_fx <- matrix( BL_boundary_adj$VarB_fx, nrow = grid_length,  ncol = grid_length, byrow = FALSE )

  ## Specify Plot Limits ##

  # If the limits for the plots aren't defined in the function call, then we assess them here by taking the range
  # over the values to be plotted.

  # Limits for the function/mean plots.
  if( identical( zlim_f, "assessed" ) ){ zlim_f <- range( c( fxu, EB_fx ) ) }
  # Limits for the variance plots.
  if( identical( zlim_var, "assessed" ) ){ zlim_var <- range( VarB_fx ) }

  ## Generate Boundary Coordinates ##

  # Boundary list for plotting on a cube using the draw_cube function.
  boundary_list <- list()

  # Use the boundary_for_plot function to generate the relevant coordinates for plotting the diagnostic slice on the cube.
  boundary_list[[1]] <- KBE::boundary_for_plot( fixed_dimension = fixed_dimension, fixed_value = fixed_value, ranges = ranges )

  # Use the boundary_for_plot function to generate the relevant coordinates for plotting the boundaries on the cube.
  boundary_list[[2]] <- KBE::boundary_for_plot( fixed_dimension = K_d, fixed_value = xK_d, ranges = ranges )
  if( identical( L_d, NA ) == FALSE ){ boundary_list[[3]] <- KBE::boundary_for_plot( fixed_dimension = L_d, fixed_value = xL_d, ranges = ranges ) }
  if( identical( M_d, NA ) == FALSE ){ boundary_list[[4]] <- KBE::boundary_for_plot( fixed_dimension = M_d, fixed_value = xM_d, ranges = ranges ) }

  ## Draw Cube ##

  # Use the draw_cube function to draw the top row cube visually depicting the known boundaries and boundary over which the subsequent plots are.
  KBE::draw_cube( coloured_hp = boundary_list,
                  lwd = lwd,
                  lty2 = lty2,
                  cex = cex,
                  col_line_width = col_line_width,
                  main = main_cube,
                  cex.main = cex.main_cube,
                  main.line = main.line )

  ## Colour Specification. ##

  # model and emulator expectation colours
  diagm_cols <- function( n ){ viridis::magma( n, begin = 0.1 ) }

  # emulator standard deviation colours
  diags_cols  <- grDevices::colorRampPalette( c( "white", "cyan", "blue", "purple", "pink" ), space = "Lab" )

  # the colours for diagnostics
  beg1 <- 0.25
  diagd_cols <- function( n ){  c( viridis::viridis( floor( n / 2 ) + 1, direc = 1, begin = beg1 )[-( floor( n / 2 ) + 1 )],
                                   viridis::plasma( ceiling( n / 2 ), direc = -1, begin = beg1 ) ) }

  ## Contour Plots.  ##

  # Split the range over which plotting is to take place into "nlevels - 1" equal-sized and "pretty" segments.
  nlevels = 20
  levels_f <- pretty( zlim_f, nlevels )
  levels_var <- pretty( zlim_var, nlevels )
  levels_diag <- pretty( zlim_diag - 0.25, nlevels ) + 0.25 # So that zero is in the centre of an interval.

  # Parameter lists for the plotting functions below.
  z_list <- list( mdxv, EB_fx, VarB_fx, ( mdxv - EB_fx ) / sqrt( VarB_fx ) )
  levels_list <- list( levels_f, levels_f, levels_var, levels_diag )
  colours_list <- list( diagm_cols, diagm_cols, diags_cols, diagd_cols )

  # Now we generate four plots.
  for( i in 1:4 ){ KBE::contour_plot( x = x1r,
                                      y = x2r,
                                      z = z_list[[i]],
                                      colours = colours_list[[i]],
                                      levels = levels_list[[i]] ) }

  ## Add Legend ##

  # Add a legend on the right hand side, if required.
  if( identical( legend, TRUE ) ){

    graphics::par( mar = c(1, 1, 1, 3.3) )

    for( i in 1:4 ){ KBE::legend_generation( colours = colours_list[[i]],
                                             levels = levels_list[[i]] ) }

  }

}
