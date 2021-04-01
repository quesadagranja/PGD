laplace_PGD <- function(n, mlim, bc, tol, maxiter, src) {
  #########################
  ###  INITIALIZATIONS  ###
  #########################
  
  # Get coordinates
  coor <- list(
    x = matrix(seq(mlim[[1]][1], mlim[[1]][2], length.out = n), nrow = 1),
    y = matrix(seq(mlim[[2]][1], mlim[[2]][2], length.out = n), nrow = 1)
  )
  # Definition of stiffness matrices
  k <- list(
    x = matrix(0, nrow = n, ncol = n),
    y = matrix(0, nrow = n, ncol = n)
  )
  # Definition of mass matrices
  m <- list(
    x = matrix(0, nrow = n, ncol = n),
    y = matrix(0, nrow = n, ncol = n)
  )
  # Definition of F matrix
  f <- list(
    x = NA,
    y = NA
  )
  # Initialization of source term matrix
  a <- list(
    x = matrix(NA, nrow = n, ncol = 0),
    y = matrix(NA, nrow = n, ncol = 0)
  )
  # Initialization of v
  v <- list()
  
  #############################
  ###  BOUNDARY CONDITIONS  ###
  #############################
  
  # Non-boundary conditions
  non_bc <- setdiff(1:n, bc)
  
  #####################################
  ###  STIFFNESS AND MASS MATRICES  ###
  #####################################
  
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
      Nx <- matrix(0, nrow = n, ncol = ngp)
      Nx[ii,] <- (gB - E2x) / gL
      Nx[ii+1,] <- (E2x - gA) / gL
      # Values of dN(x) at x = E2x for ii and ii+1 elements
      # See B in eq. 5.20 of Belytschko & Fish
      dNx <- matrix(0, nrow = n, ncol = ngp)
      dNx[ii,] <- (-1/gL) * matrix(1, nrow = 1, ncol = ngp)
      dNx[ii+1,] <- (1/gL) * matrix(1, nrow = 1, ncol = ngp)
      # Assembly and Gauss weight loop
      for (jj in 1:ngp) {
        # Mass matrix
        m[[xy]] <- m[[xy]] + gWei[jj] * (gL/2) * (Nx[,jj] %*% t(Nx[,jj]))
        # Stiffness matrix
        k[[xy]] <- k[[xy]] + gWei[jj] * (gL/2) * (dNx[,jj] %*% t(dNx[,jj]))
      }
    }
  }
  
  #####################
  ###  SOURCE TERM  ###
  #####################
  
  # Axis selection loop
  for (xy in 1:2) {
    # Evaluate the source term in all the nodes
    src_val <- src[[xy]](coor[[xy]])
    # Compose matrix a
    for (ii in 1:length(src_val)) {
      a[[xy]] <- cbind(a[[xy]], matrix(src_val[[ii]], nrow = n, ncol = 1))
    }
    # Mass matrix times source matrix
    v[[xy]] <- m[[xy]] %*% a[[xy]];
  }
  
  browser()
  
  
  # Cell memory preallocation
  
  browser()
  
  # v_size = size(source,2);
  # % Axis selection loop
  # for xy = 1:2
  # % Columns of source loop
  # for ii=1:Vsize;
  # % Evaluate the source term in all the nodes
  # A{xy}(:,ii) = source{xy,ii}(coor{xy});
  # end
  # % Mass matrix times source matrix
  # V{xy} = M{xy}*A{xy};
  # end
  
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
# # Example 1
# src <- function(x_in, y_in) {
#   x <- cos(2*pi*x_in)
#   y <- sin(2*pi*y_in)
#   return(list(x,y))
# }

# Example 3
src <- list(
  x = function(x) list(2*x^2, x, 1, 1, 3*x),
  y = function(y) list(1, 1, y^2, -0.2*y, y)
)

### CALL FUNCTION
o <- laplace_PGD(n, mlim, bc, tol, maxiter, src)