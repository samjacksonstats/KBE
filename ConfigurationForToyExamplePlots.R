### Setup for 3D Toy Example ###

## R Libraries ##

library( "viridis" )
library( "pdist" )
library( "KBE" )

## 3D Example PlotGen Function ##

# This function allows specification of a reduced argument version of Article_Plots that fixes all variables
# that will be the same across all of the plots we wish to generate.
PlotGen <- function( K_d, L_d = NA, M_d = NA,
                     xK_d, xL_d = NA, xM_d = NA ){
  
  KBE::Article_Plots( f = function( x ){ sin( x[1] / ( exp( x[2] ) ) ) + cos( x[3] ) },  # toy function
                      ranges = matrix(c( -2*pi, 2*pi, # ranges for the three parameters.
                                         -pi/4, pi/4,
                                         -2*pi, 2*pi), ncol = 2, byrow = TRUE), 
                      K_d = K_d,
                      L_d = L_d,
                      M_d = M_d,
                      xK_d = xK_d,
                      xL_d = xL_d,
                      xM_d = xM_d,
                      fixed_dimension = c( 2, 2, 1 ), # Set fixed dimensions of interest.
                      fixed_value = c( 0, -pi/8, -pi ), # Set fixed values for those variables.
                      theta = c( pi, pi/8, pi ), # Specify the correlation length parameters.
                      s2 = 2, # Specify the variance parameter.
                      zlim_f = c(-2.5, 2.5), # Specify the ranges for the function/mean...
                      zlim_var = c(0, 2), # ...and variance plots.
                      main_cube = list( bquote( x[2] == 0 ), bquote( x[2] == -pi/8 ), bquote( x[1] == -pi ) ) # set labels.
                    )

}

## 3D Example Plot Generation ##

# Single boundary: K: (x2,x3) = (0,0)
PlotGen( K_d = c(2,3), xK_d = c(0,0) )

# Two parallel known boundaries: K: (x2,x3) = (0,0), L: (x2,x3) = (0,-pi)
PlotGen( K_d = c(2,3), xK_d = c(0,0), L_d = c(2,3), xL_d = c(0,-pi) )

# Three known boundaries: K: (x2,x3) = (0,0), L: (x2,x3) = (0,-pi), M: x1 = 0
PlotGen( K_d = c(2,3), xK_d = c(0,0), L_d = c(2,3), xL_d = c(0,-pi), M_d = 1, xM_d = 0 )

# We can also consider two perpendicular boundaries: K: (x2,x3) = (0,0), L: x1 = 0
PlotGen( K_d = c(2,3), xK_d = c(0,0), L_d = 1, xL_d = 0 )

# And three perpendicular boundaries: K: x1 = 0, L: x2 = pi/8, M: x3 = 0
PlotGen( K_d = 1, xK_d = 0, L_d = 2, xL_d = pi/8, M_d = 3, xM_d = 0 )

