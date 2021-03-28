laplace_PGD <- function(n, mlim, bc, tol, maxiter, src) {
  ### INITIALIZATIONS
  # Get coordinates
  coor <- list(
    x = seq(mlim[[1]][1], mlim[[1]][2], length.out = n),
    y = seq(mlim[[2]][1], mlim[[2]][2], length.out = n)
  )
  # Definition of stiffness matrices
  k <- list(
    x = rep(0, n),
    y = rep(0, n)
  )
  # Definition of mass matrices
  m <- list(
    x = rep(0, n),
    y = rep(0, n)
  )
  # Definition of F matrix
  f <- list(
    x = NA,
    y = NA
  )
  
  ### BOUNDARY CONDITIONS
  # Non-boundary conditions
  non_bc <- setdiff(1:n, bc)
  
  ### STIFFNESS AND MASS MATRICES
  # Gauss points and weights
  gPoi <- c(-3^(-1/2), 3^(-1/2)) 
  gWei <- c(1, 1)
  # Number of Gauss points
  ngp <- length(gPoi)
  # Axis selection loop
  for (xy in 1:2) {
    # Integration loop
    for (ii in 1:(n-1)) {
      # Physical domain of integration of the current element
      gA <- coor[[xy]][ii]
      gB <- coor[[xy]][ii+1]
      gL <- gB - gA
      # Translation of gPoi from E-coordinates to x-coordinates
      E2x <- gA*(1-gPoi)/2 + gB*(gPoi+1)/2
      # Values of N(x) at x = E2x for ii and ii+1 elements
      # See N in eq. 5.20 of Belytschko & Fish
      Nx <- matrix(0, ncol = ngp, nrow = n)
      Nx[ii,] <- (gB - E2x) / gL
      Nx[ii+1,] <- (E2x - gA) / gL
      # Values of dN(x) at x = E2x for ii and ii+1 elements
      # See B in eq. 5.20 of Belytschko & Fish
      dNx <- matrix(0, ncol = ngp, nrow = n)
      dNx[ii,] <- (-1/gL) * matrix(1, ncol = ngp, nrow = 1)
      dNx[ii+1,] <- (1/gL) * matrix(1, ncol = ngp, nrow = 1)
      # Assembly and Gauss weight loop
      for (jj in 1:ngp) {
        # Mass matrix
        m[xy] = m[xy] + gWei[jj] * (gL/2) * (Nx[,jj] * T(Nx[,jj]))
        # Stiffness matrix
        k[xy] = k[xy] + gWei[jj] * (gL/2) * (dNx[,jj] * T(dNx[,jj]))
      }
    }
  }
  browser()
}
  
### INPUT PARAMETERS
# Number of nodes
n <- 11
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