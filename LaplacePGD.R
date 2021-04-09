laplace_PGD <- function(src, n, mlim, bc, tol, maxiter) {
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
    x = matrix(0, nrow = n, ncol = 0),
    y = matrix(0, nrow = n, ncol = 0)
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
    v[[xy]] <- m[[xy]] %*% a[[xy]]
  }
  
  # Size of v
  v_size <- length(src_val)
  
  ##########################################
  ###  ENRICHMENT AND PROJECTION STAGES  ###
  ##########################################
  r <- list()

  # Initialization of F error
  f_err <- 1;
  aprt <- 0;
  f_iter <- 0;

  # F loop
  while ((f_err > tol) & (f_iter < maxiter[[1]])) {
    
    ##########################
    ###  ENRICHMENT STAGE  ###
    ##########################

    # Initialization of R
    r_iter <- 0;
    r[[1]] <- matrix(rep(0, n), nrow = n, ncol = 1)
    r[[2]] <- matrix(rep(1, n), nrow = n, ncol = 1)
    r[[2]][bc] <- 0
    # Initialization of R error
    r_err <- 1
    
    # R loop
    while ((r_err > tol) & (r_iter < maxiter[[2]])) {
      # New iteration
      r_iter <- r_iter + 1
      # Copy of R
      r_aux <- r
      
      ### Get Rx when Ry is known
      
      # Left term
      d1 <-
        as.numeric(t(r[[2]]) %*% m[[2]] %*% r[[2]]) * k[[1]] + 
        as.numeric(t(r[[2]]) %*% k[[2]] %*% r[[2]]) * m[[1]]
      # Source term (V)
      d2 <- 0
      for (ii in 1:v_size) {
        d2 <- d2 + v[[1]][,ii] %*% (t(r[[2]]) %*% v[[2]][,ii])
      }
      # Right term (F)
      if (f_iter) {
        for (ii in 1:f_iter) {
          d2 <- d2 - (k[[1]] %*% f[[1]][,ii]) %*%
            (t(r[[2]]) %*% m[[2]] %*% f[[2]][,ii])
          d2 <- d2 - (m[[1]] %*% f[[1]][,ii]) %*%
            (t(r[[2]]) %*% k[[2]] %*% f[[2]][,ii])
        }
      }
      # Final assembly with imposition of domain boundaries
      r[[1]][non_bc] <- solve(d1[non_bc, non_bc]) %*% d2[non_bc]
      # Normalization
      r[[1]] <- r[[1]] / as.numeric(sqrt(t(r[[1]]) %*% m[[1]] %*% r[[1]]))
      # r[[1]] <- r[[1]] %/% norm(r[[1]])
      
      ### Get Ry when Rx is known
      
      # Left term
      d1 <- 
        as.numeric(t(r[[1]]) %*% k[[1]] %*% r[[1]]) * m[[2]] + 
        as.numeric(t(r[[1]]) %*% m[[1]] %*% r[[1]]) * k[[2]]
      # Source term (V)
      d2 <- 0
      for (ii in 1:v_size) {
        d2 <- d2 + v[[2]][,ii] %*% (t(r[[1]]) %*% v[[1]][,ii])
      }
      # Right term (F)
      if (f_iter) {
        for (ii in 1:f_iter) {
          d2 <- d2 - (m[[2]] %*% f[[2]][,ii]) %*%
            (t(r[[1]]) %*% k[[1]] %*% f[[1]][,ii])
          d2 <- d2 - (k[[2]] %*% f[[2]][,ii]) %*%
            (t(r[[1]]) %*% m[[1]] %*% f[[1]][,ii])
        }
      }
      
      # Final assembly with imposition of domain boundaries
      r[[2]][non_bc] <- solve(d1[non_bc, non_bc]) %*% d2[non_bc]
      # Normalization
      r[[2]] <- r[[2]] / as.numeric(sqrt(t(r[[2]]) %*% m[[2]] %*% r[[2]]))
      # r[[2]] <- r[[2]] %/% norm(r[[2]])
      
      # R error calculation
      r_err <- sqrt(norm(r_aux[[1]] - r[[1]]) + norm(r_aux[[2]] -  r[[2]]))
    }
    
    # New finished iteration
    f_iter <- f_iter + 1
    # Appending r to f
    for (xy in 1:2) {
      f[[xy]] <- cbind(f[[xy]], r[[xy]])
    }
    
    ##########################
    ###  PROJECTION STAGE  ###
    ##########################
    d1 <- matrix(rep(0, f_iter^2), nrow=f_iter, ncol=f_iter)
    d2 <- matrix(rep(0, f_iter), nrow=f_iter, ncol=1)
    if (f_iter) {
      for (ii in 1:f_iter) {
        # Left term (FKF, FMF)
        for (jj in 1:f_iter) {
          d1[ii, jj] <-
            (t(f[[1]][,ii]) %*% k[[1]] %*% f[[1]][,jj]) %*%
            (t(f[[2]][,ii]) %*% m[[2]] %*% f[[2]][,jj]) +
            (t(f[[1]][,ii]) %*% m[[1]] %*% f[[1]][,jj]) %*%
            (t(f[[2]][,ii]) %*% k[[2]] %*% f[[2]][,jj])
        }
        # Right term (FV)
        for (jj in 1:v_size) {
          d2[ii] <- d2[ii] + 
            (t(f[[1]][,ii]) %*% v[[1]][,jj]) %*%
            (t(f[[2]][,ii]) %*% v[[2]][,jj])
        }
      }
      
      alpha <- solve(d1) %*% d2

      for (ii in 1:f_iter) {
        sqrt_alpha <- sqrt(abs(alpha[ii]))
        f[[1]][,ii] <- sqrt_alpha * f[[1]][,ii]
        f[[2]][,ii] <- (alpha[ii] / sqrt_alpha) * f[[2]][,ii]
      }
    }
    
    # F error calculation
    f_err <-
      norm(as.matrix(f[[1]][, f_iter]), type = "2") *
      norm(as.matrix(f[[2]][, f_iter]), type = "2")
    aprt <- max(aprt, sqrt(f_err))
    f_err <- sqrt(f_err) / aprt;
  }
  
  return(f)
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

# Example 3
src <- list(
  x = function(x) list(2*x^2, x, 1, 1, 3*x),
  y = function(y) list(1, 1, y^2, -0.2*y, y)
)

### CALL FUNCTION
o <- laplace_PGD(src, n, mlim, bc, tol, maxiter)