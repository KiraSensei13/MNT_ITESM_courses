#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno Gonzalez Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     MidtermEvaluation.R
#*
#* DESCRIPTION :
#*     Simulations (Ene 19 Gpo 1)
#*     Midterm Evaluation
#*
#* NOTES :
#*     -
#*
#* START DATE :
#*     20 Mar 2019
#*******************************************************************************
################################################################################
# Second section (40 points)
# Algorithm Coding. All algorithms must be coded by yourself. The only libraries accepted for this part are the ones that already come with the R language. Do not use any package that already have the algorithm coded. Also, the code must be commented with the instructions to use your code and at least one function to demonstrate its functionality. It is recommended to do a flow chart to explain the logic of the algorithm.

#-------------------------------------------------------------------------------
#c) Runge-Kutta. The algorithm will receive an ODE, initial values for x and y, the step size and the upper bound. The output will be the plot of the second, third and fourth order RK approximations. To demonstrate the functionality of your code, use the analytical answer of the ODE to compare all the approximations (30 points)

# ----- CLEARING EVERYTHING ----

rm(list = ls(all = TRUE)) # Delete workspace
graphics.off() # Clear plots
cat("\014") # Clear console

# ----- ENTER YOUR INPUT ----

# ODE to be evaluated (enter it as a string):
funct = "-0.16*x + 0.08*x*y"
# Enter Y initial Value:
init_y = 4
# Enter X initial Value:
init_x = 4
# Enter Upper bound
upper_bound = 10
# Enter number of steps (the step size is going to be calculated calculated):
number_of_steps = 24 #length(init_x:upper_bound) * 2

# ------- THE PROGRAM ----

# Runge-Kutta - 2nd order
rungeKutta2 <- function(funct, x0, y0, x1, n) {
  f <- function(P, N) {
    return(eval(parse(text = funct), envir = list(x = P, y = N)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0) / n
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + k2) / 2
  }
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 3rd order
rungeKutta3 <- function(funct, x0, y0, x1, n) {
  f <- function(P, N) {
    return(eval(parse(text = funct), envir = list(x = P, y = N)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0) / n
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 <- h * f(x + h, y - k1 + 2 * k2)
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + 4 * k2 + k3) / 6
  }
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 4th order
rungeKutta4 <- function(funct, x0, y0, x1, n) {
  f <- function(P, N) {
    return(eval(parse(text = funct), envir = list(x = P, y = N)))
  }
  
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  h <- (x1 - x0) / n
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 <- h * f(x + 0.5 * h, y + 0.5 * k2)
    k4 <- h * f(x + h, y + k3)
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  return(cbind(vx, vy))
}

# ----- VALIDATING RESULTS ----

library (deSolve)
RKPlot <- function(func, x0, y0, x1, n) {
  # Plot of the second, third and fourth order RK approximations
  #
  # Parameters
  # ----------
  # func : string
  # ODE function of two variabless
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
  
  x2 <- rk2[, 1]
  y2 <- rk2[, 2]
  
  x3 <- rk3[, 1]
  y3 <- rk3[, 2]
  
  x4 <- rk4[, 1]
  y4 <- rk4[, 2]
  
  # Compute analytical answer of the ODE to compare all the approximations
  model <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      dy = eval(parse(text = funct), envir = list(x, y))#2*x^3 + y
      list(dy)
    })
  }
  y <- c(y = init_y)
  parms <- c()
  x <-
    seq(init_x,
        upper_bound,
        length(init_x:upper_bound) / number_of_steps)
  out <- ode(y, times = x, model, parms)
  
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
    legend = c(
      "RK 4th order",
      "RK 3rd order",
      "RK 2nd order",
      "deSolve ode function"
    ),
    col = c('#2196f3', '#4caf50', '#f44336', 'black'),
    lwd = c(1),
    bty = 'n',
    cex = .75
  )
}

# -----

# Calling the function for the solution:
RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)

#-------------------------------------------------------------------------------
#d) Coupled Runge-Kutta. Modify the fourth order RK algorithm to receive two ODEs. The output will be the plot of both ODEs approximations. (10 points)
    
        # ENTER YOUR INPUT #
# ODEs to be evaluated (enter it as a string):
funct_1 = "-0.16*x + 0.08*x*y"
funct_2 = "4.5*y - 0.9*x*y"
# Enter Y initial Value:
init_y = 4
# Enter X initial Value:
init_x = 4
# Enter Upper bound
upper_bound = 10
# Enter number of steps (the step size is going to be calculated calculated):
number_of_steps = 1000

rm(list = ls(all = TRUE)) # Delete workspace
graphics.off() # Clear plots
cat("\014") # Clear console

# -----

# Runge-Kutta - 4th order - for a 2-equation 1st order system
rungeKutta4_2eq <- function(fun_1, fun_2, x0, y0, x1, n) {
  f <- function(t, x, y) {
    return(eval(parse(text = fun_1), envir = list(
      t = t, x = x, y = y
    )))
  }
  g <- function(t, x, y) {
    return(eval(parse(text = fun_2), envir = list(
      t = t, x = x, y = y
    )))
  }
  
  xx <- double(n + 1)
  xx[1] <- xn <- x0

  yy <- double(n + 1)
  yy[1] <- yn <- y0
  
  zz <- double(n + 1)
  zz[1] <- yn <- y0

  h <- (x1 - x0) / n
  for (i in 1:n) {
    # https://www.calvin.edu/~scofield/courses/m231/materials/rungeKuttaFormulas.pdf
    tn = i
    
    kn1 = f(tn, xn, yn)
    ln1 = g(tn, xn, yn)
    
    kn2 = f(tn + h/2, xn + (1 / 2) * kn1 * h, yn + (1 / 2) * ln1 * h)
    ln2 = g(tn + h/2, xn + (1 / 2) * kn1 * h, yn + (1 / 2) * ln1 * h)
    
    kn3 = f(tn + h/2, xn + (1 / 2) * kn2 * h, yn + (1 / 2) * ln2 * h)
    ln3 = g(tn + h/2, xn + (1 / 2) * kn2 * h, yn + (1 / 2) * ln2 * h)
    
    kn4 = f(tn + h, xn + kn3 * h, yn + ln3 * h)
    ln4 = g(tn + h, xn + kn3 * h, yn + ln3 * h)
    
    xx[tn + 1] <- xn <- x0 + tn * h
    yy[tn + 1] <- xn <- xn + h * (kn1 + 2 * kn2 + 2 * kn3 + kn4)/6
    zz[tn + 1] <- yn <- yn + h * (ln1 + 2 * ln2 + 2 * ln3 + ln4)/6
    
  }
  return(cbind(xx, yy, zz))
}

out <-
  rungeKutta4_2eq("-0.16*x + 0.08*x*y", "4.5*y - 0.9*x*y", 4, 4, 24, 1000)

out <-
  rungeKutta4_2eq("-(2/25)*x + (1/50)*y", "(2/25)*x - (2/25)*y", 0, 25, 100, 1000)

plot(
  out[, 1],
  out[, 3],
  xlab = "t",
  ylab = "f(t)",
  type = "l",
  main = ' Runge-Kutta ',
  col = "black",
  lwd = 1
)

lines(out[, 1], out[, 2], col = "#4caf50")

# Calling the function for the solution:
RKPlot2(funct_1,
        funct_2,
        init_x,
        init_y,
        upper_bound,
        number_of_steps,
        -0,
        100)

################################################################################
# Third Section (60 points)
# Applied Computing. In this section you will use your fourth-order coupled Runge-Kutta algorithm code to solve different problems. Compare your answers with the "ode" function from the "deSolve" R package. Plot your answer vs the deSolve answer. Justify if the problem is stiff or not, the chosen step size, and the time units that the step size represent (seconds, minutes, hours, days, etc).

#-------------------------------------------------------------------------------
# c) Lotka-Volterra equations / predator-prey equations. The following ODEs describe the interaction between foxes and bunnies.
#
# dx/dt = -16x + 0.08xy
# dy/dt = 4.5y - 0.9xy
#
# At the start of the model there are 4 foxes and 4 bunnies. Choose the best step size to model the first 2 years of the ecosystem (30 points)

#-------------------------------------------------------------------------------
# OUR ANSWER

# ODEs to be evaluated (provided by the problem):
funct_1 = "-16*x + 0.08*x*y"
funct_2 = "4.5*y - 0.9*x*y"
# Enter Y initial Value:
init_y = 4
# Enter X initial Value:
init_x = 4
# Enter Upper bound
upper_bound = 24 # 2 yrs = 24 months
# Enter number of steps (the step size is going to be calculated calculated):
number_of_steps = 1700
 
# Calling the function for the solution:
RKPlot2(funct_1,
        funct_2,
        init_x,
        init_y,
        upper_bound,
        number_of_steps,
        -10000,
        3000)

#-------------------------------------------------------------------------------
# THE deSolve ANSWER

library(deSolve)
library(lattice)

# P -> x ; N -> y

predpreyLV <- function(t, y, p) {
  N <- y[1]
  P <- y[2]
  with(as.list(p), {
    dNdt <- c * N - d * P * N
    dPdt <- -a * P + b * P * N
    return(list(c(dNdt, dPdt)))
  })
}

a <- 0.16
b <- 0.08
c <- 4.5
d <- 0.9

p <- c(a = a,
       b = b,
       c = c,
       d = d)
y0 <- c(N = 4, P = 4)
times <- seq(0, 24, 0.1)
LV.out <- ode(y = y0, times, predpreyLV, p)

matplot(LV.out[, 1], LV.out[, 2:3], type = "l", ylab = "population")

#-------------------------------------------------------------------------------
# d) Coupled water tanks. Two water tanks are joined as depicted in the Figure 1. The model that describes the water flow from one tank to the other [gallon/minute] is given by the following ODEs:
#
# dx_1/dt = -(2/25)x_1 + (1/50)x_2
# dx_2/dt =  (2/25)x_1 - (2/25)x_2
#
# Choose the best step size to model the first 60 minutes. Assume that the first tank has 25 gallons and the second tank is empty before the h2 door is opened. (30 points)

#-------------------------------------------------------------------------------
# OUR ANSWER

# # ODEs to be evaluated (enter it as a string):
# funct_1 = "-(2/25)*x + (1/50)*y"
# funct_2 = "(2/25)*x - (2/25)*y"
# # Enter Y initial Value:
# init_y = 25
# # Enter X initial Value:
# init_x = 0
# # Enter Upper bound
# upper_bound = 60 # minutes
# # Enter number of steps (the step size is going to be calculated calculated):
# number_of_steps = 10
#
# # Calling the function for the solution:
# RKPlot2(funct_1,
#         funct_2,
#         init_x,
#         init_y,
#         upper_bound,
#         number_of_steps,
#         -150,
#         65)

#-------------------------------------------------------------------------------
# THE deSolve ANSWER

library(deSolve)
library(lattice)

# P -> x ; N -> y

watertanksLV <- function(t, y, p) {
  N <- y[1]
  P <- y[2]
  with(as.list(p), {
    dNdt <- c * P - d * N
    dPdt <- -a * P + b * N
    return(list(c(dNdt, dPdt)))
  })
}

a <- (2 / 25)
b <- (1 / 50)
c <- (2 / 25)
d <- (2 / 25)

p <- c(a = a,
       b = b,
       c = c,
       d = d)
y0 <- c(N = 0, P = 25)
times <- seq(0, 100, 0.1)
LV.out <- ode(y = y0, times, watertanksLV, p)

matplot(LV.out[, 1], LV.out[, 2:3], type = "l", ylab = "gallons")
