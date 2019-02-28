### ODEs 
## Euler method


f.x.y <- expression(log(x))
y <- expression(-x + x*log(x))
upperBound = 10
h <- 0.5
{
  
  yi <- 1
  x <- 0.001
  
  aproximaciones <- c("x" = x, "Approx" = yi)
  
  for(i in 1:(upperBound/h)){
    yi <- yi + h*eval(f.x.y)
    x <- x + h
    aproximaciones <- rbind(aproximaciones, c(x, yi))
  }
  
  print(aproximaciones)
  
  x <- seq(0,upperBound,h)
  yreal <- eval(y)
  
  
  plot(x,yreal)
  points(aproximaciones[,1],aproximaciones[,2], col = "blue", pch = 17)
}


# Heun
{
  
  yi <- 1
  x <- 0.001
  
  aproximaciones <- c("x" = x, "Approx" = yi)
  
  for(i in 1:(upperBound/h)){
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
  
  x <- seq(0,upperBound,h)
  yreal <- eval(y)
  
  
  plot(x,yreal , type = "l", col = "red")
  points(aproximaciones[,1],aproximaciones[,2], col = "blue", pch = 17)
}

#
# install.packages("deSolve")
library (deSolve)

model <- function(x, y, parms ){
  with(as.list(c(y,parms)), {
    dy = log(x)
    list(dy)
  })
}

# Initial condition
y <- c(y = 1)

parms <- c()
x <- seq(0.001,upperBound,h)
out <- ode( y, times = x, model, parms )
plot(out)

y = expression(-x+x*log(x))
x <- seq(0.001,upperBound,h)
yreal <- eval(y)

lines(x,yreal, col = "red", pch = 17)