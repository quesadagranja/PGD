library(pgd)

# ### INPUT PARAMETERS
# # Source example
# src <- list(
#   x = function(x) list(2*x^2, x, 1, 1, 3*x),
#   y = function(y) list(1, 1, y^2, -0.2*y, y)
# )
# # Number of nodes
# n <- list(
#   x = 41,
#   y = 41
# )
# # Mesh limits
# mlim <- list(
#   x = c(-1, 1),
#   y = c(-1, 1)
# )

# Source example
src <- list(
  x = function(x) list(cos(2*pi*x)),
  y = function(y) list(sin(2*pi*y))
)
# Number of nodes
n <- list(
  x = 41,
  y = 41
)
# Mesh limits
mlim <- list(
  x = c(-1, 1),
  y = c(-1, 1)
)

### CALL FUNCTION
o <- pgd::poisson_2D(src, n, mlim)

persp(o$coor$x, o$coor$y, o$t, theta=-35, phi=34, xlab="x", ylab="y", zlab="T(x,y)", ltheta=-35, lphi=55, shade=2.5, col="darkkhaki")
