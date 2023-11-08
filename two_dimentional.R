# function g to be maximized
g <- function(x,y) {
  return(-x^2 - x^2 * y^2 - 2 * x * y + 2 * x + 2)
}

# partial derivative x for function g
deriv_x <- function(x, y) {
  return(-2 * x - 2 * x * y^2 - 2 * y + 2)
}

# partial derivative y for function g
deriv_y <- function(x, y) {
  return(-2 * x^2 * y - 2 * x)
}

# gradient for function g
gradient <- function (x, y) {
  return(c(deriv_x(x,y), deriv_y(x,y)))
}

# second partial derivative x for function g
deriv_x_11 <- function(x, y) {
  return(-2 - 2 * y^2)
}

# second partial derivative x/y for function g
deriv_x_12 <- function(x, y) {
  return(-4 * x * y - 2)
}

# second partial derivative y for function g
deriv_y_22 <- function(x, y) {
  return(-2 * x^2)
}

# hessian matrix for function g
hessian_g <- function(x, y) {
  return(matrix(c(deriv_x_11(x,y), deriv_x_12(x,y), deriv_x_12(x,y), deriv_y_22(x,y)), nrow = 2, ncol = 2))
}

# produce a contour plot
xgrid <- seq(-3,3, by=0.05)
ygrid <- seq(-3,3, by=0.05)
length_x <- length(xgrid)
length_y <- length(ygrid)
dxy <- length_x * length_y
gxy <- matrix(rep(NA,dxy), nrow = length_x)
for ( i in 1:length_x) {
  for ( j in 1:length_y) {
    gxy[i, j] <- g(xgrid[i], ygrid[j])
  }
}

mgxy <- matrix(gxy, nrow = length_x, ncol = length_y)
contour(xgrid, ygrid, mgxy, nlevels=40)

# newton function
newton_g <- function(x,y, eps=0.0001) {
  xt <- c(x,y)
  xt1 <- c(x,y) + 2

  while( t(xt1 -xt) %*% (xt1 -xt) > eps) {
    xt1 <- xt
    xt <- xt1 - solve(hessian_g(xt[1] , xt[2])) %*% gradient(xt[1],xt[2])
  }
  return(xt)
}

# check definiteness of a hessian matrix
check_definiteness <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  if (all(eigenvalues > 0)) {
    cat("This Hessian Matrix is positive definite.\n")
  } else if (all(eigenvalues < 0)) {
    cat("This Hessian Matrix is negative definite.\n")
  } else {
    cat("This Hessian Matrix is neither positive nor negative definite.\n")
  }
}

# using different starting points for newton method
cat("========================================================")
converges_point1 <- newton_g(2,0)
cat("1) using starting points: (x,y) =(2,0), the newton method coverges to point: \n")
converges_point1

grad_point1 <- gradient(converges_point1[1],converges_point1[2])
cat("the corresponding gradient is: \n")
grad_point1

hessian_point1 <- hessian_g(converges_point1[1],converges_point1[2])
cat("the corresponding hessian matrix is: \n")
hessian_point1

check_definiteness(hessian_point1)

cat("========================================================")

converges_point2 <- newton_g(-1,-2)
cat("2) using starting points: (x,y) =(-1,-2) \n")
converges_point2

grad_point2 <- gradient(converges_point2[1],converges_point2[2])
cat("the corresponding gradient is: \n")
grad_point2

hessian_point2 <- hessian_g(converges_point2[1],converges_point2[2])
cat("the corresponding hessian matrix is: \n")
hessian_point2

check_definiteness(hessian_point2)

cat("========================================================")

converges_point3 <- newton_g(0,1)
cat("3) using starting points: (x,y) =(0,1) \n")
converges_point3

grad_point3 <- gradient(converges_point3[1],converges_point3[2])
cat("the corresponding gradient is: \n")
grad_point3

hessian_point3 <- hessian_g(converges_point3[1],converges_point3[2])
cat("the corresponding hessian matrix is: \n")
hessian_point3

check_definiteness(hessian_point3)

cat("========================================================")

converges_point4 <- newton_g(1,-1)
cat("4) using starting points: (x,y) =(1,-1) \n")
converges_point4

grad_point4 <- gradient(converges_point4[1],converges_point4[2])
cat("the corresponding gradient is: \n")
grad_point4

hessian_point4 <- hessian_g(converges_point4[1],converges_point4[2])
cat("the corresponding hessian matrix is: \n")
hessian_point4

check_definiteness(hessian_point4)





