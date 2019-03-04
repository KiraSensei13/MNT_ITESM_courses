# library("cubature")

##################
#Exercise 1:
##################

library(rSymPy)

sympy("var('x')")
sympy("(sin(x)).series(x, 0, 10)")
# Fórmula: sum( x^(2*k+1)/ factorial(2*k + 1))

sympy("(exp(I*x)).series(x, 0, 10)")
# Fórmula: cos(x) + I*sin(x)


##################
# Exercise 2 and 3:
##################
# Lower and Upper bound
a <- -2
b <- 2

# Taylor expression
sympy("(exp(x + x**2)).series(x, 0, 6)")

# Riemann sum
deltaX <- 0.01
y <- expression(exp(x + x^2))

# Generate the steps. 
# NOTE: Remember that we do not take the last value of the upper bound
x <- seq(a,b-deltaX,deltaX)
# Evaluation of y at each step
y.x <- eval(y)
# Multiply by deltaX
y.x.delta <- y.x * deltaX
# Sum all terms
RiemmanSum <- sum(y.x.delta)

# I print the riemann sum
RiemmanSum

# Lets compare the results
Real <- integrate( function(x){exp(x + x^2)}, a, b)

error <- abs(Real$value - RiemmanSum)

##################
# Exercise 4:
##################

z <- expression(exp(((x + y)^2)))
a <- 0
b <- 1

deltaX <- 0.01
deltaY <- 0.01

# Generate the steps. 
# NOTE: Remember that we do not take the last value of the upper bound
x.all <- seq(a,b-deltaX,deltaX)
y.all <- seq(a,b-deltaY,deltaY)
# Evaluation of z at each step
z.x.y <- matrix(0, ncol = 100, nrow = 100 )
for(i in 1:100){
  y = y.all[i]
  for(j in 1:100){
    x = x.all[j]
    z.x.y[j,i] <- eval(z)    
  }
}


# Multiply by deltaX
z.deltas <- z.x.y * deltaX * deltaY
# Sum all terms
RiemmanSum <- sum(z.deltas)

# I print the riemann sum
RiemmanSum
# Real

library(pracma)
real <- integral2( function(x,y){exp((x + y)^2)}, a, b, a, b)
# 4.89

error <- RiemmanSum - real$Q
error

##################
## Exercise 5 (Extra)
##################

delta <- 0.1
xi <- 3
y <- expression(x^3)

x <- xi - delta
y.left <-  eval(y)

x <- xi + delta
y.right <-  eval(y)

x <- xi
y.center <- eval(y)

y.dx.dx <- (y.right - 2*y.center + y.left)/(delta*delta)


# Using the real values (derivatives)
y.dx.a <- 3*(x^2) # First derivative
y.dx.dx.a <- 2*(3*(x)) # Second derivative

print(y.dx.dx)   # Approximation
print(y.dx.dx.a) # Real value


