biharmonic <- function(dataint,
                       grid_res = 100) {

  x_min <- 0
  x_max <- 1
  y_min <- 0
  y_max <- 1

  xo <- seq(x_min,
            x_max,
            length = grid_res)
  yo <- seq(y_min,
            y_max,
            length = grid_res)
  xo <- matrix(rep(xo,
                   grid_res),
               nrow = grid_res,
               ncol = grid_res)
  yo <- t(matrix(rep(yo, grid_res),
                 nrow = grid_res,
                 ncol = grid_res))

  xy_coords <- unique(dataint[, c("x", "y")])
  xy <- xy_coords[, 1, drop = TRUE] + xy_coords[, 2, drop = TRUE] * sqrt(as.complex(-1))
  d <- matrix(rep(xy,
                  length(xy)),
              nrow = length(xy),
              ncol = length(xy))
  d <- abs(d - t(d))
  diag(d) <- 1
  g <- (d^2) * (log(d) - 1) #Green's function
  diag(g) <- 0
  weights <- qr.solve(g, dataint$signal)
  xy <- t(xy)

  # Remind me to make this code readable at some point.
  outmat <- numeric(grid_res^2)

  interp_man <- function(xo, yo) {

    tmp <-
      matrix(complex(real = xo,
                     imaginary = yo),
             nrow = grid_res,
             ncol = grid_res)
    tmp_x <- numeric(length(xy))

    for (i in seq(length(tmp))) {
      tmp_x <- (abs(tmp[i] - xy)^2) * (log(abs(tmp[i] - xy)) - 1)
      tmp_x <- ifelse(is.nan(tmp_x),
                      0,
                      tmp_x)
      outmat[i] <- tmp_x %*% weights
    }
    outmat
  }

  outmat <- interp_man(xo, yo)
  data2 = data.frame(x = rep(seq(1,grid_res),times=grid_res),
                     y = rep(seq(1,grid_res),each=grid_res),
                     outmat)



  return(data2)
}
