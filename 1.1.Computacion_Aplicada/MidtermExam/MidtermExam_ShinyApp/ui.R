# Decided to develop an App in order to make the MidTerm exam more user Freindly.
# Install the shiny package:
# install.package("shiny")
# Load Library:
# library(shiny)

shinyUI(
# Defining the user interface:
  fluidPage(
    
# Adding the Header:
    
    titlePanel(h3("Midterm Exam", align = "center")),
    titlePanel(h6("by: Antonio Osamu Katagiri Tanaka (A01212611) and Bruno Gonzalez Soria (A01169284)", align = "center")),

# This is the side panel where the variables can be defined by the user to be evaluated by the program:  
  
    sidebarPanel(
      selectInput("Exam Problems", "Please Select the problem to be evaluated:",
                  choices=c("First Section - Problem 1 a)", "Second Section - Problem a)", "Second Section - Problem b)", "Second Section - Problem c)")),
      
      conditionalPanel(condition = "input.Exam Problems == 'First Section - Problem 1 a)'"),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem a)'"),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem b)'",
                       textInput("point", "Please enter the point x to which you want to approximate: ", 0),
                       textInput("xmin", "Please enter lower value of x: ", -6),
                       textInput("xmax", "Please enter upper value of x:", 6),
                       sliderInput("factors", "Please select number of factors:", min = 1, max = 5000, value = 1, step = 1)),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem c)'")
      
      
    ),
    
    mainPanel(
      
      tabsetPanel(type="tab",
                  tabPanel("Solution",
                           tags$b(paste("Instructions:")),
                           textOutput("Problem Instructions"),plotOutput("Graphical Representation")),
                  tabPanel("Code", verbatimTextOutput("Code")),
                  tabPanel("Chocolate Data", tableOutput("Chocolate"))
                  )
      
    )
    
  )
  
  
)
