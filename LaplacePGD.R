laplace_PGD <- function(n, mlim, bc, tol, maxiter, src) {
  ### INITIALIZATIONS
  # Preallocation
  K <- list(x=NA, y=NA)
  M <- list(x=NA, y=NA)
  F <- list(x=NA, y=NA)
  coor <- list(x=NA, y=NA)
  # Axis selection loop
  for (xy in 1:2) {
    # Get coordinates
    browser()
    coor[xy] <- seq(mlim[[xy]][1], mlim[[xy]][2], n)
  }
  
  # %% INITIALIZATIONS
  # % Cell memory preallocation   
  # K = cell(2,1);
  # M = cell(2,1);
  # F = cell(2,1);
  # coor = cell(2,1);
  # % Axis selection loop
  # for xy = 1:2
  # % Get coordinates
  # coor{xy} = linspace(Mlim{xy}(1), Mlim{xy}(2), n);
  # % Definition of stiffness matrices
  # K{xy} = zeros(n);
  # % Definition of mass matrices
  # M{xy} = zeros(n);
  # % Definition of F matrix
  # F{xy} = [];
  # end
}
  
### INPUT PARAMETERS
# Number of nodes
n <- 10
# Mesh limits
mlim <- list(
  x = c(-1, 1),
  y = c(-1, 1)
)
# Boundary conditions
bc <- c(1, n)
# Error tolerance
tol <- 1e-4
# Maximum iterations
maxiter <- list(
  f_loop = 500,
  r_loop = 251
)

### DEFINITION OF THE SOURCE TERM
# Example 1
src <- function(x_in, y_in) {
  x = cos(2*pi*x_in)
  y = sin(2*pi*y_in)
  return(list(x,y))
}

### CALL FUNCTION
o <- laplace_PGD(n, mlim, bc, tol, maxiter, src)