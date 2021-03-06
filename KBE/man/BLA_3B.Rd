% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BLA_3B.r
\name{BLA_3B}
\alias{BLA_3B}
\title{Bayes Linear Adjustment by 3 Known Boundaries}
\usage{
BLA_3B(
  x,
  K_d,
  L_d,
  M_d,
  xK = NA,
  xL = NA,
  xM = NA,
  xLK = NA,
  xMK = NA,
  xML = NA,
  xMLK = NA,
  xK_d = NA,
  xL_d = NA,
  xM_d = NA,
  fxK,
  fxL,
  fxM,
  fxLK,
  fxMK,
  fxML,
  fxMLK,
  E_fx = 0,
  E_fxK = 0,
  E_fxL = 0,
  E_fxM = 0,
  E_fxLK = 0,
  E_fxMK = 0,
  E_fxML = 0,
  E_fxMLK = 0,
  theta,
  s2
)
}
\arguments{
\item{x}{points at which we want to update}

\item{K_d}{the dimensions which, when fixed at certain values, result in known boundary K.}

\item{L_d}{the dimensions which, when fixed at certain values, result in known boundary L.}

\item{M_d}{the dimensions which, when fixed at certain values, result in known boundary M.}

\item{xK}{the projection of x onto known boundary K}

\item{xL}{the projection of x onto known boundary L}

\item{xM}{the projection of x onto known boundary M}

\item{xLK}{the projection of x first onto known boundary L and then known boundary K.}

\item{xMK}{the projection of x onto the intersection of known boundaries K and M.}

\item{xML}{the projection of x onto the intersection of known boundaries L and M.}

\item{xMLK}{the projection of x onto the intersection of M and that obtained by projecting first
onto L and then onto K.}

\item{xK_d}{values the dimensions K must take for the function to be known}

\item{xL_d}{values the dimensions L must take for the function to be known}

\item{xM_d}{values the dimensions M must take for the function to be known}

\item{fxK}{function evaluated at x projected onto the boundary K.}

\item{fxL}{function evaluated at x projected onto the boundary L.}

\item{fxM}{function evaluated at x projected onto the boundary M.}

\item{fxLK}{function evaluated at xLK.}

\item{fxMK}{function evaluated at xMK.}

\item{fxML}{function evaluated at xML.}

\item{fxMLK}{function evaluated at xMLK.}

\item{E_fx}{prior expectation for the function f(x)}

\item{E_fxK}{prior expectation for f(x^K)}

\item{E_fxL}{prior expectation for f(x^L)}

\item{E_fxM}{prior expectation for f(x^M)}

\item{E_fxLK}{prior expectation for f(x^LK)}

\item{E_fxMK}{prior expectation for f(x^MK)}

\item{E_fxML}{prior expectation for f(x^ML)}

\item{E_fxMLK}{prior expectation for f(x^MLK)}

\item{theta}{vector of correlation length parameter values.}

\item{s2}{scalar variance parameter value.}
}
\value{
\item{EB_fx}{Expected value of f(x) adjusted by knowledge of function behaviour along K, L and M.}
\item{VarB_fx}{Variance of f(x) adjusted by knowledge of function behaviour along K, L and M.}
\item{CovB_fx}{Covariance of f(x) adjusted by knowledge of function behaviour along K, L and M.}
}
\description{
Perform a Bayes linear adjustment utilising knowledge of function behaviour along
three known boundaries in the input space.  In this case boundaries K and L should be parallel to each other,
and M should be perpendicular to K and L.
}
\examples{
# Toy function
f <- function( x ){

sin( x[1] / ( exp( x[2] ) ) ) + cos( x[3] )

}
x <- matrix( runif( 12 ), ncol = 3 )
K_d = 2
L_d = c(2,3)
M_d = 1
xK_d = 0
xL_d = c(1,1)
xM_d = 0
#' # If we are in a parallel setting, then xLK (projection of x first onto L and then K)
# is given as follows:
xLK_d <- xL_d
xLK_d[1:length(xK_d)] <- xK_d
# And LK_d (fixed values of coordinates for projections first onto L and then K)
# is just given by L_d.
LK_d <- L_d
fxK <- f_boundary( x = x, K_d = K_d, xK_d = xK_d, f = f )
fxL <- f_boundary( x = x, K_d = L_d, xK_d = xL_d, f = f )
fxM <- f_boundary( x = x, K_d = M_d, xK_d = xM_d, f = f )
fxLK <- f_boundary( x = x, K_d = LK_d, xK_d = xLK_d, f = f )
fxMK <- f_boundary( x = x, K_d = c(K_d, M_d), xK_d = c(xK_d, xM_d), f = f )
fxML <- f_boundary( x = x, K_d = c(L_d, M_d), xK_d = c(xL_d, xM_d), f = f )
fxMLK <- f_boundary( x = x, K_d = c(LK_d, M_d), xK_d = c(xLK_d, xM_d), f = f )
theta <- c( pi, pi/8, pi )
s2 <- 2
BA <- BLA_3B( x = x, K_d = K_d, L_d = L_d, M_d = M_d,
              xK_d = xK_d, xL_d = xL_d, xM_d = xM_d,
              fxK = fxK, fxL = fxL, fxM = fxM,
              fxLK = fxLK, fxMK = fxMK, fxML = fxML, fxMLK = fxMLK,
              theta = theta, s2 = s2 )
}
