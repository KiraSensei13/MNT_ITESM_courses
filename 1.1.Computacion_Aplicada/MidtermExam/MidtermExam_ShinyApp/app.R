#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     MidtermExam.R
#*
#* DESCRIPTION :
#*     Exam (Ene 19 Gpo 1)
#*     Midterm Evaluation
#* 
#* NOTES : 
#*     -
#* 
#* START DATE : 
#*     07 Mar 2019  
#*******************************************************************************

########### INITIAL COMMENTS AND PREPARATIONS ##################################

# Decided to develop an App in order to make the MidTerm exam more user Freindly.
# Install the shiny package:
# install.packages("shiny")
# Load Library:
library(shiny)
# Giving some style:
# install.packages("shinythemes")
# Load Library:
library(shinythemes)
# Adding the codes:
library(markdown)
library(deSolve)
library(gsubfn)


#*******************************************************************************

########### USER INTERFACE #####################################################

# Defining the user interface:
ui <- fluidPage(
  
########### STYLE ############################################################## 
  
  # Defining style:
  theme = shinytheme("cerulean"),
  
########### HEADER OF THE APPLICATION ##########################################
  
  titlePanel(h1("Midterm Exam", align = "center")),
  titlePanel(h5("by: Antonio Osamu Katagiri Tanaka (A01212611) and Bruno Gonzalez Soria (A01169284)", align = "center")),
  
########### THIS IS THE SIDE PANEL FOR THE USER TO ENTER INPUT #################  
  
  sidebarPanel(
    
    ######### Selection of Problem #############################################
    
    selectInput("ExamProblems", "Please Select the problem to be evaluated:",
                choices = c("First Section - Problem 1 a)", "Second Section - Problem a)", "Second Section - Problem b)", "Second Section - Problem c)")),
    
    ######### Input Panel for Problem 1 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'First Section - Problem 1 a)'"),
    
    ######### Input Panel for Problem 2 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'Second Section - Problem a)'",
      textInput("valsX", "Enter the values for the 'x' column separated only by a comma:", "1,2,3,4,7,8,9,10"),
      textInput("valsY", "Enter the values for the 'y' column separated only by a comma:", "4,3,2,1,3,6,8,10"),
      helpText("Make sure you add the same amount of values in both columns!")
      ),
    
    ######### Input Panel for Problem 3 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'Second Section - Problem b)'",
      textInput("funT","Please enter your function:", "sin(x^2)"),
      numericInput("pnt", "Please enter the 'x' value to which you want to approximate: ", 0, step = 0.1),
      numericInput("trms", "Please enter the desired ammount of terms:", 2, min = 1 , max = 8),
      checkboxInput("show", "Show all approximations?", TRUE),
      numericInput("xmin", "Edit 'x' axis lower limit:", -6),
      numericInput("xmax", "Edit 'x' axis upper limit:", 6)
      ),
    
    ######### Input Panel for Problem 4 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'Second Section - Problem c)'",
      textInput("ODE", "Please enter the ODE", "2*x^3 + y"),
      numericInput("xi","Enter x initial value:",0),
      numericInput("yi", "Enter y initial value:",0),
      numericInput("step","Enter the desired step size:",1),
      numericInput("upbound","Enter the upper bound:",2),
      selectInput("RK", "Select which RK order you want to plot:", choices = c("RK2","RK3","RK4","All together"))
      )
    
  ),
  
########### THIS IS WHERE THE INFORMATION WILL BE SHOWN ########################
  mainPanel(
    
    # Displaying the information in tabs:
    tabsetPanel(type="tab",
                tabPanel("Solution",
                         tags$b(paste("Instructions:")),
                         textOutput("Problem Instructions"),
                         plotOutput("Graphical Representation")),
                tabPanel("Code",
                         uiOutput("Procedure")),
                tabPanel("Chocolate",
                         tableOutput("Chocolate"))
    )
  )
)

#*******************************************************************************


#*******************************************************************************

########### SERVER #############################################################

# Defining the servre back-end function:

server <- function (input,output,session){
  
########### THIS SECTION IS JUST FOR THE INSTRUCTIONS OF EACH PROBLEM ##########
  
  output$`Problem Instructions` <- renderText({
    
    problemId <- input$`ExamProblems`
    
    if (problemId == "First Section - Problem 1 a)"){
      
      instructions <- "Database: Chocolate.
      Dependent variable: Chocolate -> (Chocolate is a categorical variable)
      Choose only the most significant variables to model the dependent variable. (10 points)"
      
    }
    
    else{
      
      if (problemId == "Second Section - Problem a)"){
        
        instructions <- "Lagrange polynomials. This algorithm receives a nx2 matrix, where the first column represents the x coordinate while the second column represents the y coordinate. The code must provide as output the Lagrange polynomial interpolation expression in terms of x. (30 points) *TIP: Use the functions: expression, D, parse and paste within a loop to get the desired output."
        
      }
      
      else{
        
        if (problemId == "Second Section - Problem b)"){
          
          instructions <- "Taylor series: The algorithm receives an expression or a string with the function to do and the number of terms to get. The output will be an expression containing all the Taylor series about 0. (30 points) *TIP: Use the functions: expression, D, parse and paste within a loop to get the desired output."
          
        }
        
        else{
          
          if (problemId == "Second Section - Problem c)"){
            
            instructions <- "Runge-Kutta. The algorithm will receive an ODE, initial values for x and y, the step size and the upper bound. The output will be the plot of the second, third and fourth order RK approximations. To demonstrate the functionality of your code, use the analytical answer of the ODE to compare all the approximations (30 points)"
            
          }
          
          else {
            
            instructions <- "Sorry for the inconvenience, we are still working hard to make this application evenn better."
            
          }
        }
      }
    }
    
    paste(instructions)
    
  })
  
  
########### THIS SECTION IS TO ADD THE CODE (PROCEDURE) ########################
  
  output$`Procedure` <- renderText({
    
    problemId <- input$`ExamProblems`
    
    ################## Procedure for problem 1 ########################
    
    if (problemId == "First Section - Problem 1 a)") {
      
      proc <- includeMarkdown("Problem1.Rmd")
      
    }
    
    else{
    
    ################## Procedure for problem 2 ########################
    
      if (problemId == "Second Section - Problem a)") {
        
        proc <- includeMarkdown("Problem2.Rmd")
      
      }
    
      else{
      
    ################## Procedure for problem 3 ########################    
    
        if (problemId == "Second Section - Problem b)") {
          
          proc <- includeMarkdown("Problem3.Rmd")
      
        }
    
        else{
    
    ################## Procedure for Problem 4 ########################
    
          if (problemId == "Second Section - Problem c)"){
            
            proc <- includeMarkdown("Problem4.Rmd")
      
          }
    
          else{
            
            proc <- "Sorry for the inconvenience, we are still working hard to make this application evenn better."
      
          }
        }
      }
    }
    
    paste(proc)
    
  })
  
  
########### THIS SECTION IS FOR THE PLOTS CREATED BY THE PROGRAM ###############
  
  output$`Graphical Representation` <- renderPlot({
    
    problemId <- input$`ExamProblems`
    
    ############## Plot for Problem 1 Chocolate ################################
    
    if (problemId == "First Section - Problem 1 a)") {
      
      choco_data = read.csv("Chocolate.csv")
      chocolate_labels <- choco_data[,1] 
      choco_data <- choco_data[,-1] 
      choco_data$chocolate <- as.factor(choco_data$chocolate)
      mdl <- glm(chocolate ~ ., data = choco_data, family = "binomial")
      mdl <- glm(chocolate ~ fruity + winpercent, data = choco_data, family = "binomial")
      predicted <- predict(mdl, choco_data[,], type = "response")
      glm_predicted = ifelse(predicted > 0.5, 1, 0)
      
      graph <- hist(glm_predicted, main = "Histogram of Chocolate Frequency in Candy", ylab = "Frequency of Trait", xlab = "Presence of Chocolate", ylim = c(0, 60), labels = TRUE)
      
    }
    
    else{
      
      ############## Plot for Problem 2 Lagrange ###############################
      
      if (problemId == "Second Section - Problem a)") {
        
        library(rSymPy)
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
          
          plot(xx, yy, xlab="x", ylab="f(x)")
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
        
        x = as.numeric(unlist(strapply(strsplit(input$`valsX`,","), "\\d+", as.numeric)))
        y = as.numeric(unlist(strapply(strsplit(input$`valsY`,","), "\\d+", as.numeric)))
        graph <- lagrange(cbind(x,y))
        
      }
      
      else{
        
        ############## Plot for problem 3 Taylor ###############################
        
        if (problemId == "Second Section - Problem b)") {
          
          pnt <- input$`pnt`
          xmin <- input$`xmin`
          xmax <- input$`xmax`
          trms <- seq(1:input$`trms`)
          
          library(pracma)
          taylorPlot <- function(funct, taylorOrder) {
            # Plot the Taylor approximations up to the 2nd, 4th, 6th and 8th terms
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
            x <- seq(xmin, xmax, length.out = 250)
            yf <- f(x)
            c <- pnt # The Taylor Series shall be centered in zero
            
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
              xlim = c(xmin, xmax),
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
          
          if (input$`show` == TRUE){
            for (i in trms) {
              graph <- taylorPlot(input$`funT`, i)
            }
          }
          else{
            graph <- taylorPlot(input$`funT`, input$`trms`)
          }
          
        }
        
        else{
          
          ############## Plot for problem 4 Runge Kutta ########################
          
          if (problemId == "Second Section - Problem c)") {
            
            if (input$`RK` == "RK2" | input$`RK` == "All together"){
              
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
              
              graph <- RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)
              
            }
            else{
              
              if(input$`RK` == "RK3" | input$`RK` == "All together"){
                
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
                
                graph <- RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)
                
              }
              else{
                if(input$`RK` == "RK4" | input$`RK` == "All together"){
                  
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
                  
                  graph <- RKPlot(funct, init_x, init_y, upper_bound, number_of_steps)
                  
                }
                else{
                  paste("Sorry for the inconvenience, we are still working hard to make this application evenn better.")
                }
              }
            }
            
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
            
            funct           = input$`ODE`
            init_y          = input$`yi`
            init_x          = input$`xi`
            upper_bound     = input$`upbound`
            number_of_steps = input$`step`
            
          }
          
          else{
            
            graph <- "Sorry for the inconvenience, we are still working hard to make this application evenn better."
            
            
          }
          
        }
      }
    }
    
    paste(graph)
    
  })
  
  ########### THIS SECTION CORRESPONDS ONLY TO THE CHOCOLATE DATA ##############  
  
  output$`Chocolate` <- renderTable({
    
    choco_data = read.csv("Chocolate.csv")
    
    choco_data
    
  })
  
}

#*******************************************************************************

########### RUN THE APP ########################################################

shinyApp(ui = ui, server = server)

#*******************************************************************************
#                                                                              #
#                        THANKS FOR YOUR SUPPORT!                              #
#                                                                              #
#*******************************************************************************