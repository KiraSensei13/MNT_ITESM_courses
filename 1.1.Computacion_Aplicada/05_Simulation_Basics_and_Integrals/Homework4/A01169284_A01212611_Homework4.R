#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     Homework4.R
#*
#* DESCRIPTION :
#*     Simulations (Ene 19 Gpo 1)
#*     Homework 4
#*
#* NOTES :
#*     - https://www.math.ubc.ca/~pwalls/math-python/integration/riemann-sums/
#*
#* START DATE :
#*     21 Feb 2019
#*******************************************************************************

# Plotting the Areas under Curves ##############################################
integralPlot <- function(f, a, b, from = a, to = b, title = NULL) {
  x <- seq(from, to, length.out = 100) # input continuum
  y <- f(x) # output
  
  # plot the curve
  plot(x, y,
       xlim = c(from, to),
       ylim = c(ifelse(min(y) < 0, min(y), 0), max(y)),
       xlab = "x",
       ylab = "f(x)",
       main = title,
       col.main = "#86B875",
       type = "l",
       lwd = 3,
       col = "#86B875")
  
  # area under the curve
  x <- seq(a, b, length.out = 100)
  y <- f(x)
  polygon(c(x, b, a, a), c(y, 0, 0, f(a)),
          border = adjustcolor("#7DB0DD", alpha.f = 0.3),
          col = adjustcolor("#7DB0DD", alpha.f = 0.3))
  }

# PART 2 #######################################################################
# RIEMANN SUMS FUNCTION

riemann_sum <- function(f,a,b,n) {
  # initialize values
  lower.sum <- 0;
  upper.sum <- 0;
  h <- (b - a)/n;
  
  # riemann right sum
  for (i in n:1) {
    x <- a + i*h;
    lower.sum <- lower.sum + f(x);
  }
  lower.sum <- h*lower.sum;
  
  # riemann left sum
  for (i in 1:n) {
    x <- b - i*h;
    upper.sum <- upper.sum + f(x);
  }
  upper.sum <- h*upper.sum;
  
  # print/get riemann sum
  cat(sprintf("The true value is between %f and %f.\n",
          as.double(lower.sum),
          as.double(upper.sum)));
  return(c(lower.sum,upper.sum));
}

# -----

# let's generate some functions to test our algorithm
f0 <- function(x) {
  res = sin(x);
  return(res);
}

f1 <- function(x) {
  res = x;
  return(res);
}

f2 <- function(x) {
  res = 4/(1+x^2);
  return(res);
}

# -----

riemann_sum(f0,0,pi/2,10)
riemann_sum(f1,0,1,10000)
riemann_sum(f2,0,1,10000)


# PART 3 #######################################################################
# Integrate the function f(x)=exp(x+x^2) from -2 to 2, using Rieman sums.
f3 <- function(x) {
  res = exp(x+x^2);
  return(res);
}

# compute riemann_sum for f3
riemann_sum(f3,-2,2,100000)
# let's verify our calualtion using R's function
integrate(function(x){exp(x + x^2)}, lower=-2, upper=2)
# let's plot the curve
integralPlot(f=f3, a=-2, b=2, title=expression(exp(x+x^2)))

# plot(labx=(-2,2), function(x){exp(x+x^2)}, type="h")
# step = 0.001
# f.x = expression (exp(x + x^2))
# x = 1.01+step
# eval(f.x)*step
