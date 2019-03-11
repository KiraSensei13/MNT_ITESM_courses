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

chocolate_labels <- choco_data[,1]
choco_data <- choco_data[,-1]

# Setting numeric values to categoric value:
choco_data$chocolate <- as.factor(choco_data$chocolate)
summary(choco_data)

# Since the dependent variable "Chocolate" is a categorical variable, we have decided to use Logistic Regression:
mdl <- glm(chocolate ~ ., data = choco_data, family = "binomial")
summary(mdl)

mdl <- glm(chocolate ~ fruity + caramel +  peanutyalmondy + nougat + crispedricewafer + hard + bar + pluribus + sugarpercent + pricepercent + winpercent, data = choco_data, family = "binomial")
summary(mdl) # Visualizing the model to determine the most significant variables.

# Lets see the performance of our beta coefficients:
predicted <- predict(mdl, choco_data[,], type = "response")
hist(predicted)

#################################################################################
# SECOND SECTION
# a) Lagrange polynomials. This algorithm receives a nx2 matrix, where the first column represents the x coordinate while the second column represents the y coordinate. The code must provide as output the Lagrange polynomial interpolation expression in terms of "x". (30 points)


################################################################################
# b) Taylor series: The algorithm receives an expression or a string with the function to do and the number of terms to get. The output will be an expression containing all the Taylor series about 0. (30 points)



################################################################################
# c) Runge-Kutta. The algorithm will receive an ODE, initial values for x and y, the step size and the upper bound. The output will be the plot of the second, third and fourth order RK approximations. To demonstrate the functionality of your code, use the analytical answer of the ODE to compare all the approximations (30 points)
