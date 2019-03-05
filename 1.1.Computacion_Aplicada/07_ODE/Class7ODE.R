### ODEs 
## Euler method


f.x.y <- expression((-2*(x^3)) + (12*(x^2)) - (20*x) + 8.5)


y <- expression((-0.5*(x^4)) + (4*(x^3)) - (10*(x^2)) + (8.5*x) + 1 )



h <- 0.25
h <- 0.05
{

yi <- 1
x <- 0

aproximaciones <- c("x" = x, "Approx" = yi)

for(i in 1:(4/h)){
  
  yi <- yi + h*eval(f.x.y)
  x <- x + h
  aproximaciones <- rbind(aproximaciones, c(x, yi))
}

print(aproximaciones)

x <- seq(0,4,0.05)
yreal <- eval(y)


plot(x,yreal , 
     ylim = c(min(c(aproximaciones[,2],yreal)), 
              max(c(aproximaciones[,2],yreal))), 
     type = "l", col = "red")


points(seq(0,4,h),aproximaciones[,2], col = "blue", pch = 17)
}


#


f.x.y <- expression((-2*(x^3)) + (12*(x^2)) - (20*x) + 8.5)


y <- expression((-0.5*(x^4)) + (4*(x^3)) - (10*(x^2)) + (8.5*x) + 1 )


h <- 0.5
h <- 0.25
{

yi <- 1
x <- 0

aproximaciones <- c("x" = x, "Approx" = yi)

for(i in 1:(4/h)){
  yreal <- eval(y)
  
  y1 <- eval(f.x.y)
  y.temp <- yi 
  
  
  yi <- yi + h*y1
  
  
  
  x <- x + h
  y2 <- eval(f.x.y)
  
  
  yi <- y.temp  + (y1 + y2)*h*0.5
  aproximaciones <- rbind(aproximaciones, c(x, yi))
}

print(aproximaciones)

x <- seq(0,4,0.01)
yreal <- eval(y)


plot(x,yreal , 
     ylim = c(min(c(aproximaciones[,2],yreal)), 
              max(c(aproximaciones[,2],yreal))), 
     type = "l", col = "red")


points(seq(0,4,h),aproximaciones[,2], col = "blue", pch = 17)
}

#
# install.packages("deSolve")
library (deSolve)
#


model <- function(x, y, parms ){
  with(as.list(c(y,parms)), {
    dy = (-2*(x^3)) + (12*(x^2)) - (20*x) + 8.5
    list(dy)
  })
}

# Initial condition
y <- c(y = 1)

parms <- c()

x <- seq(0,4,0.5)


out <- ode( y, times = x, model, parms )
plot(out)


y <- expression((-0.5*(x^4)) + (4*(x^3)) - (10*(x^2)) + (8.5*x) + 1 )
x <- seq(0,4,0.01)
yreal <- eval(y)

lines(seq(0,4,0.01),yreal, col = "red", pch = 17)




# Equations:
# y'1 = ((I2 - I3)/I1)*y2*y3
# y'2 = ((I3 - I1)/I2)*y1*y3
# y'3 = ((I1 - I2)/I3)*y2*y1

model3 <- function(t,y,parms){
  with(as.list(c(y,parms)), {
    dy1 = ((I2 - I3)/I1)*y[2]*y[3]
    dy2 = ((I3 - I1)/I2)*y[1]*y[3]
    dy3 = ((I1 - I2)/I3)*y[2]*y[1]
    list(c(dy1,dy2,dy3))
  })
}

# Initial conditions:
# y1(0) = 1, y2(0) = 0, y3(0) = 0.9
yini <- c(y1 = 1, y2 = 0, y3 = 0.9)
# Constants:
parms <- c(I1 = 0.5, I2 = 2, I3 = 3)
# Independent variable
times <- seq(0,20,0.01)

# ODE
out3 <- ode(yini, times, model3, parms)
head(out3)

## Symbolic derivatives

#

exprs <- expression(x^3)

# "x"
dy.dx <- deriv(exprs, "x")
dy.dx


x <- 0:10


eval(dy.dx)

# Multivariate
f.x.y <- expression(x^2 + y^2)

gradient.f.x.y <- deriv(f.x.y, c("x","y"))
gradient.f.x.y


# Other way to do it
D(f.x.y, "x" )
D(f.x.y, "y" )


x <- 2
y <- 1
eval(gradient.f.x.y)


deriv( f.x.y, c("x","y"), func = TRUE)


exprs <- expression(x^3)
d.dx <- D(exprs, "x")
d2.dx2 <- deriv(d.dx, "x")


exprs <- expression(x^3)
d.dx <- D(exprs, "x")
d2.dx2 <- D(d.dx, "x")
d3.dx3 <- deriv(d2.dx2, "x")



f.x.y <- expression(x^2 + y^2)
d2.f.x.y <- deriv(f.x.y, c("x", "y"), hessian = T)






