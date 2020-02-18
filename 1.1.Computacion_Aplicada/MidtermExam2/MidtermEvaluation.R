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

# ----- CLEARING EVERYTHING -----

rm(list = ls(all = TRUE)) # Delete workspace
graphics.off() # Clear plots
cat("\014") # Clear console

# ----- CODE IMPLEMENTATION -----

# Runge-Kutta - 2nd order
rungeKutta2 <- function(funct, x0, y0, x1, n) {
  
  # convert the string-function into an actual function
  f <- function(xx, yy) {
    return(eval(parse(text = funct), envir = list(x = xx, y = yy)))
  }
  
  # create vx and vy to store the RK calculations for each step
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  
  # get the step size h
  h <- (x1 - x0) / n
  
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    
    # compute the RK constants k
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    
    # store the calculated x and y within vx and vy respectively
    # and update x and y
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + k2) / 2
  }
  
  # return vx and vy (just in case)
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 3rd order
rungeKutta3 <- function(funct, x0, y0, x1, n) {
  
  # convert the string-function into an actual function
  f <- function(xx, yy) {
    return(eval(parse(text = funct), envir = list(x = xx, y = yy)))
  }
  
  # create vx and vy to store the RK calculations for each step
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  
  # get the step size h
  h <- (x1 - x0) / n
  
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    
    # compute the RK constants k
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 <- h * f(x + h, y - k1 + 2 * k2)
    
    # store the calculated x and y within vx and vy respectively
    # and update x and y
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + 4 * k2 + k3) / 6
  }
  
  # return vx and vy (just in case)
  return(cbind(vx, vy))
}

# -----

# Runge-Kutta - 4th order
rungeKutta4 <- function(funct, x0, y0, x1, n) {
  
  # convert the string-function into an actual function
  f <- function(xx, yy) {
    return(eval(parse(text = funct), envir = list(x = xx, y = yy)))
  }
  
  # create vx and vy to store the RK calculations for each step
  vx <- double(n + 1)
  vy <- double(n + 1)
  vx[1] <- x <- x0
  vy[1] <- y <- y0
  
  # get the step size h
  h <- (x1 - x0) / n
  
  for (i in 1:n) {
    # from http://campus.usal.es/~mpg/Personales/PersonalMAGL/Docencia/MetNumTema4Teo(09-10).pdf
    
    # compute the RK constants k
    k1 <- h * f(x, y)
    k2 <- h * f(x + 0.5 * h, y + 0.5 * k1)
    k3 <- h * f(x + 0.5 * h, y + 0.5 * k2)
    k4 <- h * f(x + h, y + k3)
    
    # store the calculated x and y within vx and vy respectively
    # and update x and y
    vx[i + 1] <- x <- x0 + i * h
    vy[i + 1] <- y <- y + (k1 + 2 * k2 + 2 * k3 + k4) / 6
  }
  
  # return vx and vy (just in case)
  return(cbind(vx, vy))
}

# -----

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
  
  # Calculate RK aproxminations
  rk2 = rungeKutta2(func, x0, y0, x1, n) # RK 2nd order
  rk3 = rungeKutta3(func, x0, y0, x1, n) # RK 3rd order
  rk4 = rungeKutta4(func, x0, y0, x1, n) # RK 4th order
  
  # store the RK approxamations in separate variables
  x2 <- rk2[, 1] # RK 2nd order
  y2 <- rk2[, 2]
  
  x3 <- rk3[, 1] # RK 3rd order
  y3 <- rk3[, 2]
  
  x4 <- rk4[, 1] # RK 4th order
  y4 <- rk4[, 2]
  
  # Compute analytical answer of the ODE using deSolve's ode fuction
  # to compare all RK the approximations
  model <- function(x, y, parms) {
    with(as.list(c(y, parms)), {
      dy = eval(parse(text = funct), envir = list(x, y))
      list(dy)
    })
  }
  y <- c(y = init_y)
  parms <- c()
  x <-
    seq(init_x,
        upper_bound,
        length(init_x:upper_bound) / number_of_steps)
  # save deSolve's results
  out <- ode(y, times = x, model, parms)
  
  # Plot deSolve's output
  plot(
    out,
    xlab = "x",
    ylab = "f'(x)",
    type = "l",
    main = ' Runge-Kutta ',
    col = "black",
    lwd = 2
  )
  # Plot RK approximations
  lines(x2, y2, col = "#f44336") # RK 2nd order (red)
  lines(x3, y3, col = "#4caf50") # RK 3rd order (green)
  lines(x4, y4, col = "#2196f3") # RK 4th order (blue)
  
  # add a legent to the plot
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

# ----- ENTER YOUR INPUT -----

# ODE to be evaluated:
funct = "x^2 - y"
# Enter Y initial Value:
init_y = 1
# Enter X initial Value:
init_x = -5
# Enter Upper bound
upper_bound = 5
# Enter number of steps (the step size is going to be calculated calculated):
number_of_steps = length(init_x:upper_bound) * 2

# ----- VALIDATING RESULTS -----

# Calling the function for the solution:
RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)


#-------------------------------------------------------------------------------
#d) Coupled Runge-Kutta. Modify the fourth order RK algorithm to receive two ODEs. The output will be the plot of both ODEs approximations. (10 points)

# ----- CLEARING EVERYTHING -----

rm(list = ls(all = TRUE)) # Delete workspace
graphics.off() # Clear plots
cat("\014") # Clear console

# ----- CODE IMPLEMENTATION -----

# Runge-Kutta - 4th order - for a 2-equation 1st order system
rungeKutta4_2eq <- function(y0, times, funct, p) {
  
  # create xx and yy to store the RK calculations for each step
  xx <- double(length(times) - 1)
  xx[1] <- xn <- y0[1]
  
  yy <- double(length(times) - 1)
  yy[1] <- yn <- y0[2]
  
  # get the step size h
  h <- (max(times) - min(times)) / length(times)
  i = 2 # start populating the 2nd element of xx and yy. the 1st ones are y0
  
  for (tn in times) {
    # from https://www.calvin.edu/~scofield/courses/m231/materials/rungeKuttaFormulas.pdf
    # and https://math.stackexchange.com/questions/721076/help-with-using-the-runge-kutta-4th-order-method-on-a-system-of-2-first-order-od
    
    # compute the RK constants k and l
    kn1 = funct(tn, xn, yn)[[1]][[1]]
    ln1 = funct(tn, xn, yn)[[1]][[2]]
    
    kn2 = funct(
      tn + h / 2,
      xn + (1 / 2) * kn1 * h,
      yn + (1 / 2) * ln1 * h)[[1]][[1]]
    ln2 = funct(
      tn + h / 2,
      xn + (1 / 2) * kn1 * h,
      yn + (1 / 2) * ln1 * h)[[1]][[2]]
    
    kn3 = funct(
      tn + h / 2,
      xn + (1 / 2) * kn2 * h,
      yn + (1 / 2) * ln2 * h)[[1]][[1]]
    ln3 = funct(
      tn + h / 2,
      xn + (1 / 2) * kn2 * h,
      yn + (1 / 2) * ln2 * h)[[1]][[2]]
    
    kn4 = funct(tn + h, xn + kn3 * h, yn + ln3 * h)[[1]][[1]]
    ln4 = funct(tn + h, xn + kn3 * h, yn + ln3 * h)[[1]][[2]]
    
    # store the calculated xn and yn within xx and yy respectively
    xx[i] = as.numeric(xn)
    yy[i] = as.numeric(yn)
    
    # update xn and yn
    xn <- xn + h * (kn1 + 2 * kn2 + 2 * kn3 + kn4) / 6
    yn <- yn + h * (ln1 + 2 * ln2 + 2 * ln3 + ln4) / 6
    
    # update the index i to avoid overwriting in xx and yy
    i = i + 1
  }
  
  # remove the last element within xx and yy,
  # their length shall match times's length to plot
  xx = xx[-length(times + 1)]
  yy = yy[-length(times + 1)]
  
  # mmin and mmax are the plot y-limits
  mmin = as.numeric(min(xx))
  mmax = as.numeric(max(xx))
  mmin <- if (mmin > min(yy))
    min(yy)
  else
    mmin
  mmax <- if (mmax < max(yy))
    max(yy)
  else
    mmax
  
  # plot xx
  plot(
    times,
    xx,
    xlab = "t",
    ylab = "f(t)",
    type = "l",
    main = ' Runge-Kutta ',
    col = "#2196f3", # (blue)
    lwd = 1,
    ylim = c(mmin, mmax)
    
  )
  # plot yy in the same image
  lines(times, yy, col = "#4caf50") # (green)
  
  # add a legent to the plot
  legend(
    'topleft',
    inset = .05,
    legend = c(
      "f1",
      "f2"
    ),
    col = c('#4caf50','#2196f3'),
    lwd = c(1),
    bty = 'n',
    cex = .75
  )
  
  # return xx and yy (just in case)
  return(cbind(xx, yy))
}

# ----- ENTER YOUR INPUT -----

# ODEs to be evaluated: P -> x ; N -> y ;
myFunction <- function(t, y, p) {
  N <- y[1]
  P <- y[2]
  with(as.list(p), {
    dNdt <- c * N - d * P * N
    dPdt <- -a * P + b * P * N
    return(list(c(dNdt, dPdt)))
  })
}

# ODEs' coefficients:
a <- 0.16
b <- 0.08
c <- 4.5
d <- 0.9

# ODEs' coefficients/parameters:
p <- c(a = a,
       b = b,
       c = c,
       d = d)

# Initial conditions: P -> x ; N -> y ;
y0 <- c(N = 4, P = 4)

# time vector: to define t0, tf, and the step size
times <- seq(0, 12, 0.1) # time unit = months

# ----- VALIDATING RESULTS -----

out <- rungeKutta4_2eq(y0, times, myFunction, p)


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

# P -> x ; N -> y

# ----- ENTER YOUR INPUT -----

# ODEs to be evaluated: P -> x ; N -> y ;
predpreyLV <- function(t, y, p) {
  N <- y[1]
  P <- y[2]
  with(as.list(p), {
    dNdt <- c * N - d * P * N
    dPdt <- -a * P + b * P * N
    return(list(c(dNdt, dPdt)))
  })
}

# ODEs' coefficients:
a <- 16
b <- 0.08
c <- 4.5
d <- 0.9

# ODEs' coefficients/parameters:
p <- c(a = a,
       b = b,
       c = c,
       d = d)

# Initial conditions: P -> x ; N -> y ;
y0 <- c(N = 4, P = 4)

# time vector: to define t0, tf, and the step size
times <- seq(0, 2, 0.01) # time unit = years

# ----- STUDENT RESULTS -----

out <- rungeKutta4_2eq(y0, times, predpreyLV, p)

# ----- deSolve RESULTS -----

library(deSolve)
library(lattice) # to plot using matplot

LV.out <- ode(y = y0, times, predpreyLV, p)
matplot(LV.out[, 1], LV.out[, 2:3], type = "l",main = 'deSolve Approximation', xlab = "t", ylab = "population")

# ----- RESULTS ANALYSIS -----

# This problem is stiff as the differential equations are not numerically stable and must use excesively small steplenghts due to the rapid changes in population that force the entire system to be evaluated with such small steplengths, even in the regions where the curve is smoother. We decided to use the stepsize of 0.01 years (87.66 hours), since this stepsize proved to have a stable response even with larger time inetvals.


#-------------------------------------------------------------------------------
# d) Coupled water tanks. Two water tanks are joined as depicted in the Figure 1. The model that describes the water flow from one tank to the other [gallon/minute] is given by the following ODEs:
#
# dx_1/dt = -(2/25)x_1 + (1/50)x_2
# dx_2/dt =  (2/25)x_1 - (2/25)x_2
#
# Choose the best step size to model the first 60 minutes. Assume that the first tank has 25 gallons and the second tank is empty before the h2 door is opened. (30 points)

# P -> x ; N -> y

# ----- ENTER YOUR INPUT -----

# ODEs to be evaluated: P -> x ; N -> y ;
watertanksLV <- function(t, y, p) {
  N <- y[1]
  P <- y[2]
  with(as.list(p), {
    dNdt <- c * P - d * N
    dPdt <- -a * P + b * N
    return(list(c(dNdt, dPdt)))
  })
}

# ODEs' coefficients:
a <- (2 / 25)
b <- (1 / 50)
c <- (2 / 25)
d <- (2 / 25)

# ODEs' coefficients/parameters:
p <- c(a = a,
       b = b,
       c = c,
       d = d)

# Initial conditions: P -> x ; N -> y ;
y0 <- c(N = 0, P = 25)

# time vector: to define t0, tf, and the step size
times <- seq(0, 600, 1) # time unit = minutes

# ----- STUDENT RESULTS -----

out <- rungeKutta4_2eq(y0, times, watertanksLV, p)

# ----- deSolve RESULTS -----

library(deSolve)
library(lattice) # to plot using matplot

LV.out <- ode(y = y0, times, watertanksLV, p)
matplot(LV.out[, 1], LV.out[, 2:3], type = "l", main = 'deSolve Approximation', xlab = "t", ylab = "gallons")

# ----- RESULTS ANALYSIS -----

# This is not a stiff problem, since the numerical method applied (coupled-RK) has a stable behaviour, well appoximated to the solution and doesn't need to be evaluated with excesively small steplengths in relation to the smoothness of the curve. The behaviour continues to be stable even in larger time intervals. We decided to use the steplength of 0.1 minutes (6 seconds) since it is a comprehensive length and to avoid the bit of noise generated at the beginning of the function when using biger steps. Nonetheless, increasing the time interval to be evaluated allows an increase in steplength with equaly good performance.

#*******************************************************************************

