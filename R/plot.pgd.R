plot.pgd <- function(o, nmodes = dim(o$f$x)[2], phi = 45, theta = 45) {
  # Compute surface
  z <- o$f$x[,1:nmodes] %*% t(o$f$y[,1:nmodes])
  # Plot surface
  persp(
    x = o$coor$x,
    y = o$coor$y,
    z = z,
    phi = phi,
    theta = theta,
    xlab = "x",
    ylab = "y"
  )
}