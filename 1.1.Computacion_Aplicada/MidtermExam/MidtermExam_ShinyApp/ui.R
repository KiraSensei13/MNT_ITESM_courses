shinyUI(
  pageWithSidebar(
    headerPanel(title = "Midterm Exam"),
    
    sidebarPanel(
      selectInput("Exam Problems", "Please Select the problem to be evaluated:",
                  choices=c("First Section - Problem 1 a)", "Second Section - Problem a)", "Second Section - Problem b)", "Second Section - Problem c)")),
      
      conditionalPanel(condition = "input.Exam Problems == 'First Section - Problem 1 a)'"),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem a)'"),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem b)'"),
      
      conditionalPanel(condition = "input.Exam Problems == 'Second Section - Problem c)'")
      
      
    ),
    
    mainPanel(
      
      tags$b(paste("Instructions:")),
      textOutput("Problem Instructions"),
      plotOutput("Graphical Representation")
      
    )
    
  )
  
  
)
