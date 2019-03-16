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
      textInput("valsX", "Enter the values for the 'x' column separated by a space:", "1 2 3 4"),
      textInput("valsY", "Enter the values for the 'y' column separated by a space:", "4 3 2 1")
      ),
    
    ######### Input Panel for Problem 3 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'Second Section - Problem b)'",
      textInput("funT","Please enter your function:", "sin(x^2)"),
      numericInput("pnt", "Please enter the 'x' value to which you want to approximate: ", 0),
      numericInput("xmin", "Please enter lower value of 'x': ", -6),
      numericInput("xmax", "Please enter upper value of 'x':", 6),
      numericInput("trms", "Please enter the desired ammount of terms:", 2),
      checkboxInput("show", "Show all approximations?", TRUE)
      ),
    
    ######### Input Panel for Problem 4 ########################################
    
    conditionalPanel(
      condition = "input.ExamProblems == 'Second Section - Problem c)'",
      textInput("ODE", "Please enter the ODE", "x+y"),
      numericInput("xi","Enter x initial value:",0),
      numericInput("yi", "Enter y initial value:",0),
      numericInput("step","Enter the desired step size:",1),
      numericInput("upbound","Enter the upper bound:",2),
      selectInput("RK", "Select which RK order you want to plot:", choices = c("RK1","RK2","RK3","RK4","All together"))
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
    
    ############## Plot for Problem 1 #################################
    
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
      
      ############## Plot for Problem 2 ###############################
      
      if (problemId == "Second Section - Problem a)") {
        
        
        
      }
      
      else{
        
        ############## Plot for problem 3 #############################
        
        if (problemId == "Second Section - Problem b)") {
          
          library(pracma)
          pnt <- input$`pnt`
          xmin <- input$`xmin`
          xmax <- input$`xmax`
          trms <- seq(1,input$`trms`,1)
          fctn <- input$`funT`
          npoints <- 250
          nmult <- 5
          
          #####
          
          poly <- function(a,p) {
            
            A <- matrix(0,4,5)
            A[1,1:2] <- c(a[1]-a[2]*p,a[2])
            A[2,1:3] <- A[1,1:3]+a[3]*c(p^2,-2*p,1)
            A[3,1:4] <- A[2,1:4]+a[4]*c(p^3,-3*p^2,3*p^2,1)   
            A[4,] <- A[3,]+a[5]*c(p^4,-4*p^3,6*p^2,-4*p,1)
            A
          }
          
          data <- reactive({
            x <- seq(xmin,xmax,length=npoints)
            if(pnt==0) i <- nmult*25
            else i <- nmult*pnt
            p <- x[i]
            h <- x[2]-x[1]
            f <- function(x) {
              eval(parse(text=fctn))
            }
            y <- f(x)
            p0 <- y[i]
            p1 <- (y[i+1]-y[i])/h
            p2 <- (y[i-1]-2*y[i]+y[i+1])/h^2
            p3 <- (y[i+2]-3*y[i+1]+3*y[i]-y[i-1])/h^3
            p4 <- (y[i+2]-4*y[i+1]+6*y[i]-4*y[i-1]+y[i-2])/h^4
            
            yr <- c(min(y)-(max(y)-min(y)/3),max(y)+(max(y)-min(y)/3))
            
            list(cbind(x,y),c(p0,p1,p2,p3,p4),p,yr)
          })
          
          x <- data()[[1]][,1]
          y <- data()[[1]][,2]
          yr <- data()[[4]]
          
          graph <-  plot(x,y,ylim=yr, xlim = c(xmin,xmax), xlab="x",ylab="",type="l",lwd=3)
                    if(pnt==0){
                      i <- nmult*25
                    } 
                    else{
                      i <- nmult*pnt
                    } 
                    points(x[i],y[i],pch=20,cex=2)
                    a <- data()[[2]]
                    p <- data()[[3]]
                    
                    y <- a[1]+a[2]*(x-p)
                    lines(x,y,lwd=1,col="blue")
                     
                    y <- y+a[3]/2*(x-p)^2
                    lines(x,y,lwd=1,col="green")
                     
                    y <- y+a[4]/6*(x-p)^3
                    lines(x,y,lwd=1,col="red")
                     
                    y <- y+a[5]/24*(x-p)^4
                    lines(x,y,lwd=1,col="gray")
                    
        }
        
        else{
          
          ############## Plot for problem 4 ###########################
          
          if (problemId == "Second Section - Problem c)") {
            
            ODE <- input$`ODE`
            xi <- input$`xi`
            yi <- input$`yi`
            stp <- input$`step`
            upbnd <- input$`upbound`
            RKn <- input$`RK`
            
            
            
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