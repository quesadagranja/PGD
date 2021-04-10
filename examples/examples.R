library(pgd)

### INPUT PARAMETERS
# Number of nodes
n <- list(
  x = 41,
  y = 41
)
# Boundary conditions
bc <- list(
  x = c(1, n$x),
  y = c(1, n$y)
)
# Mesh limits
mlim <- list(
  x = c(-1, 1),
  y = c(-1, 1)
)
# Error tolerance
tol <- 1e-3
# Maximum iterations
maxiter <- list(
  f_loop = 500,
  r_loop = 251
)

# Source example
src <- list(
  x = function(x) list(2*x^2, x, 1, 1, 3*x),
  y = function(y) list(1, 1, y^2, -0.2*y, y)
)

### CALL FUNCTION
o <- pgd::laplacian_2D(src, n, bc, mlim, tol, maxiter)

### PLOT FUNCTION
plot(o)