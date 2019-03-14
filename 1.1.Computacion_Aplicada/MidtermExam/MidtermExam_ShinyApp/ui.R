######## INITIAL COMMENTS AND PREPARATIONS #####################################

# Decided to develop an App in order to make the MidTerm exam more user Freindly.
# Install the shiny package:
# install.package("shiny")
# Load Library:
# library(shiny)
# Giving some style:
# install.packages("shinythemes")
# Load Library:
# library(shinythemes)

shinyUI(
  
######## STYLE##################################################################
# Defining the user interface:
  fluidPage(
    
# Defining style:
    theme = shinytheme("sandstone"),
    
######## HEADER OF THE APPLICATION #############################################
    
    titlePanel(h1("Midterm Exam", align = "center")),
    titlePanel(h5("by: Antonio Osamu Katagiri Tanaka (A01212611) and Bruno Gonzalez Soria (A01169284)", align = "center")),

######## THIS IS THE SIDE PANEL FOR THE USER TO ENTER INPUT ####################  
  
    sidebarPanel(
      selectInput("ExamProblems", "Please Select the problem to be evaluated:",
                  choices = c("First Section - Problem 1 a)", "Second Section - Problem a)", "Second Section - Problem b)", "Second Section - Problem c)")),
      
      conditionalPanel(
        condition = "input.ExamProblems == 'First Section - Problem 1 a)'"),
      
      conditionalPanel(
        condition = "input.ExamProblems == 'Second Section - Problem a)'"),
      
      conditionalPanel(
        condition = "input.ExamProblems == 'Second Section - Problem b)'",
        numericInput("point", "Please enter the point x to which you want to approximate: ", 0),
        numericInput("xmin", "Please enter lower value of x: ", -6),
        numericInput("xmax", "Please enter upper value of x:", 6),
        sliderInput("factors", "Please select number of factors:",
                    min = 1, max = 5000, value = 1, step = 1)),
      
      conditionalPanel(
        condition = "input.ExamProblems == 'Second Section - Problem c)'")
      
      
    ),
    
####### THIS IS WHERE THE INFORMATION WILL BE SHOWN ############################
    mainPanel(

# Displaying the information in tabs:
      tabsetPanel(type="tab",
                  tabPanel("Solution",
                           tags$b(paste("Instructions:")),
                           textOutput("Problem Instructions"),
                           plotOutput("Graphical Representation")),
                  tabPanel("Code", verbatimTextOutput("Code")),
                  tabPanel("Chocolate",
                           tableOutput("Chocolate"))
                  )
      
              )
      
    )
  
  
)
