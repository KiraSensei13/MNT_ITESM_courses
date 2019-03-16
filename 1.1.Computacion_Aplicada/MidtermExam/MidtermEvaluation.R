#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     MidtermExam.R
#*
#* DESCRIPTION :
#*     Simulations (Ene 19 Gpo 1)
#*     Midterm Evaluation
#*
#* NOTES :
#*     - 
#*
#* START DATE :
#*     07 Mar 2019
#*******************************************************************************
################################################################################
# FIRST SECTION
# Database: Chocolate
# a) Dependent variable: Chocolate -> Chocolate is a categorical variable
# Choose only the most significant variables to model the dependent variable. (10 points)

# First we define a variable to call for the database to be analized:
choco_data = read.csv("Chocolate.csv")
str(choco_data) # Verifying there is no repeated information.
summary(choco_data) # Exploring the data.

# There is no need to filter our data since all variables are categorical, with the exceptions of "sugarpercent", "pricepercent" and "winpercent" but none of them have any row with a value of 0. The only information that causes irregularities is the competitors name.

chocolate_labels <- choco_data[,1] # First column contains the names of candy and we define it as the label for each row.
choco_data <- choco_data[,-1] # Taking away the first column that contains the candy names.

# Setting numeric values to categoric value:
choco_data$chocolate <- as.factor(choco_data$chocolate)
summary(choco_data)

# Since the dependent variable "Chocolate" is a categorical variable, we have decided to use Logistic Regression:
mdl <- glm(chocolate ~ ., data = choco_data, family = "binomial")
summary(mdl)

mdl <- glm(chocolate ~ fruity + caramel +  peanutyalmondy + nougat + crispedricewafer + hard + bar + pluribus + sugarpercent + pricepercent + winpercent, data = choco_data, family = "binomial")
summary(mdl) # Visualizing the model to determine the most significant variables.

# Takng away all the not significative parameters we obtain the following model:

mdl <- glm(chocolate ~ fruity + winpercent, data = choco_data, family = "binomial")
summary(mdl) # Visualizing the model to determine the most significant variables.

# Lets see the performance of our beta coefficients:
predicted <- predict(mdl, choco_data[,], type = "response")
hist(predicted)

#We assume that: all values less than 0.5 correspond to the first category, 0.5 or bigger correspond to the second category
glm_predicted = ifelse(predicted > 0.5, 1, 0)
hist(glm_predicted)

################################################################################
# SECOND SECTION
# a) Lagrange polynomials. This algorithm receives a nx2 matrix, where the first column represents the x coordinate while the second column represents the y coordinate. The code must provide as output the Lagrange polynomial interpolation expression in terms of "x". (30 points) *TIP: Use the functions: expression, D, parse and paste within a loop to get the desired output.

library(polynom)
library(pracma)
lagrange <- function(coordinates) {
  # Plot the Lagrange polynomial, evaluated within the x coordinates
  # (from parameter coordinates)
  #
  # Parameters
  # ----------
  # coordinates : nx2 matrix 
  # where the first column represents the x coordinate while the second
  # column represents the y coordinate
  #
  # Returns
  # -------
  # The Lagrange polynomial

  x = coordinates[,1]
  y = coordinates[,2]
  
  interPoly = poly.calc(x,y)
  xx = x
  yy = lagrangeInterp(x,y,xx)
  
  plot(xx, yy, xlab="x", ylab="f(x)", main = ' Lagrange polynomials ')
  lines(interPoly, col = "#4caf50")
  
  legend(
    'topleft',
    inset = .05,
    legend = c("Coordinates", "Lagrange polynomial"),
    col = c('black', '#4caf50'),
    lwd = c(1),
    bty = 'n',
    cex = .75
  )
  
  return(interPoly)
}

# -----

x = seq(-5,5,0.5)
y = x^3 + 5*x^2 + 1
lagrange(cbind(x,y))

################################################################################
# b) Taylor series: The algorithm receives an expression or a string with the function to do and the number of terms to get. The output will be an expression containing all the Taylor series about 0. (30 points) *TIP: Use the functions: expression, D, parse and paste within a loop to get the desired output.

library(pracma)
taylorPlot <- function(funct, taylorOrder) {
  # Plot the Taylor approximations up to the taylorOrder terms
  #
  # Parameters
  # ----------
  # f : function
  # Vectorized function of one variable
  # taylorOrder : numeric
  # the number of terms to get
  #
  # Returns
  # -------
  # The Taylor Series
  
  f <- function(n) {
    return(eval(parse(text=funct), envir=list(x=n)))
  }
  
  # Interval of points to be ploted
  from = -2*pi
  to = 2*pi
  x <- seq(from, to, length.out = 100)
  yf <- f(x)
  c <- 0 # The Taylor Series shall be centered in zero
  
  taylorOut <- taylor(f, c, taylorOrder)
  yp <- polyval(taylorOut, x)
  
  plot(
    x,
    yf,
    xlab = "x",
    ylab = "f(x)",
    type = "l",
    main = ' Taylor Series Approximation of f(x) ',
    col = "black",
    lwd = 2
  )
  
  lines(x, yp, col = "#4caf50")
  
  legend(
    'topleft',
    inset = .05,
    legend = c("Taylor Approximation", "f(x)"),
    col = c('#4caf50', 'black'),
    lwd = c(1),
    bty = 'n',
    cex = .75
  )
  
  return(taylorOut)
}

# -----

taylorPlot("sin(x)", 8)
taylorPlot("cos(x)", 2)
taylorPlot("tan(x)", 3)

################################################################################
# c) Runge-Kutta. The algorithm will receive an ODE, initial values for x and y, the step size and the upper bound. The output will be the plot of the second, third and fourth order RK approximations. To demonstrate the functionality of your code, use the analytical answer of the ODE to compare all the approximations (30 points)

# Runge-Kutta - 2nd order
rungeKutta2 <- function(funct, x0, y0, x1, n) {
  f <- function(xx,yy) {
    return(eval(parse(text=funct), envir=list(x=xx,y=yy)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0)/n
  for(i in 1:n) {
    k1 <- h*f(x, y)
    k2 <- h*f(x + 0.5*h, y + 0.5*k1)
    vx[i + 1] <- x <- x0 + i*h
    vy[i + 1] <- y <- y + (k1 + k2)/2
  }
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 3rd order
rungeKutta3 <- function(funct, x0, y0, x1, n) {
  f <- function(xx,yy) {
    return(eval(parse(text=funct), envir=list(x=xx,y=yy)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0)/n
  for(i in 1:n) {
    k1 <- h*f(x, y)
    k2 <- h*f(x + 0.5*h, y + 0.5*k1)
    k3 <- h*f(x + 0.5*h, y + 0.5*k2)
    vx[i + 1] <- x <- x0 + i*h
    vy[i + 1] <- y <- y + (k1 + 4*k2 + k3)/6
  }
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 4th order
rungeKutta4 <- function(funct, x0, y0, x1, n) {
  f <- function(xx,yy) {
    return(eval(parse(text=funct), envir=list(x=xx,y=yy)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0)/n
  for(i in 1:n) {
    k1 <- h*f(x, y)
    k2 <- h*f(x + 0.5*h, y + 0.5*k1)
    k3 <- h*f(x + 0.5*h, y + 0.5*k2)
    k4 <- h*f(x + h, y + k3)
    vx[i + 1] <- x <- x0 + i*h
    vy[i + 1] <- y <- y + (k1 + 2*k2 + 2*k3 + k4)/6
  }
  return(cbind(vx, vy))
}

# -----

RKPlot <- function(func, x0, y0, x1, n) {
  # Plot of the second, third and fourth order RK approximations
  #
  # Parameters
  # ----------
  # func : string
  # function of two variable
  # x0, y0 : numeric
  # initial values
  # x1 : numeric
  # upper bound
  # n : numeric
  # number of steps
  #
  # Returns
  # -------
  # void
  
  # Calculate aproxminations
  rk2 = rungeKutta2(func, x0, y0, x1, n)
  rk3 = rungeKutta3(func, x0, y0, x1, n)
  rk4 = rungeKutta4(func, x0, y0, x1, n)
  
  x2 <- rk2[,1]
  y2 <- rk2[,2]
  
  x3 <- rk3[,1]
  y3 <- rk3[,2]
  
  x4 <- rk4[,1]
  y4 <- rk4[,2]
  
  # Computeanalytical answer of the ODE to compare all the approximations
  model <- function(x, y, parms){
    with(as.list(c(y,parms)), {
      dy = eval(parse(text=funct), envir=list(x,y))#2*x^3 + y
      list(dy)
    })
  }
  y <- c(y = init_y)
  parms <- c()
  x <- seq(init_x,upper_bound,length(init_x:upper_bound)/number_of_steps)
  out <- ode( y, times = x, model, parms )
  
  # Plot
  plot(
    out,
    xlab = "x",
    ylab = "f'(x)",
    type = "l",
    main = ' Runge-Kutta ',
    col = "black",
    lwd = 2
    )

  lines(x2, y2, col = "#f44336")
  lines(x3, y3, col = "#4caf50")
  lines(x4, y4, col = "#2196f3")
  
  legend(
    'topleft',
    inset = .05,
    legend = c("RK 4th order", "RK 3rd order", "RK 2nd order", "f'(x)"),
    col = c('#2196f3', '#4caf50', '#f44336', 'black'),
    lwd = c(1),
    bty = 'n',
    cex = .75
  )
}

# -----

funct           = "x^2 - 0.05*y"
init_y          = 1
init_x          = -10
upper_bound     = 10
number_of_steps = length(init_x:upper_bound)*2

RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)

