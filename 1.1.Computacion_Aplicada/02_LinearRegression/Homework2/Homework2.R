#************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     Homework2.R
#*
#* DESCRIPTION :
#*     Computación aplicada (Ene 19 Gpo 1)
#*     Homework 2
#*
#* NOTES :
#*     - https://www.zoology.ubc.ca/~schluter/R/fit-model/
#*
#* START DATE :
#*     23 Jan 2019
#************************************************************************

#install.packages('e1071', dependencies=TRUE)
#install.packages('caret', dependencies=TRUE)

#*** PART 1 *************************************************************

# Do a regression analysis to the “pima” dataset from the “faraway” library
library(faraway)
library(caret)
data(pima)

# * Analyze the database and select only the observations with no missing
#   data
# complete.cases() returns a logical vector with the value TRUE for rows that are complete, and FALSE for rows that have some NA values
completeData = complete.cases(pima)
# remove rows with incomplete data
mydata = pima[completeData, ]

#Let's take a look to the data ...
str(mydata)
pairs(mydata)
cor(mydata)
summary(mydata)

# * Use the R function to fit a model.
# * Dependent variable : Test -> test is a categorical variable ...
# Create the regression model.
#names(mydata)

#LOGISTIC REGRESSION
glmFitModel <- glm(test ~ pregnant + glucose + diastolic + triceps + insulin + bmi + diabetes + age, data = mydata)
summary(glmFitModel)
#Let's remove triceps, insulin, and age as they do not have a significant impact on the model
glmFitModel <- glm(test ~ pregnant + glucose + diastolic + bmi + diabetes, data = mydata)
summary(glmFitModel)
glmFitModel

#get test predictions using the LOGISTIC REGRESSION
glm_predictions = predict(glmFitModel, mydata)
#hist(glm_predictions)

#We assume that: all values less than 0.5 correspond to the first category, 0.5 or bigger correspond to the second category
glm_predictions = ifelse(glm_predictions > 0.5, 1, 0)

# get the percentage of prediction ...
variableTable = table(glm_predictions , as.character(mydata$test))
sum(diag(variableTable))/sum(variableTable)

#count how many of the predictions are right vs wrong
correct_lm_predictions = mydata$test == glm_predictions
#FALSE:WRONG TRUE:RIGHT
summary(correct_lm_predictions)

# * Confusion matrix to see the model performance
confusionMatrix(as.factor(glm_predictions), as.factor(mydata$test), positive = NULL, dnn = c("Prediction", "Reference"))

# Which regression did you use?
# * ANCOVA, ANOVA, simple regression, logistic regression
# * Justify your answer
sprintf("Logistic regression, since logistic regression is better used when the dependent variable is categorical.")


#*** PART 2 *************************************************************

# Do a regression analysis to the “teengamb” dataset from the “faraway”
# library
library(faraway)
data(teengamb)
# complete.cases() returns a logical vector with the value TRUE for rows
# that are complete, and FALSE for rows that have some NA values
completeData = complete.cases(teengamb)
# remove rows with incomplete data
mydata = teengamb[completeData, ]

#Let's take a look to the data ...
str(mydata)
pairs(mydata)
cor(mydata)
summary(mydata)

#There is data from 2 sex, labeled 0, 1
#Let's convert the sex data into a factor
mydata$sex = factor(as.numeric(mydata$sex))

# * Use the normal equations to fit the 𝛽 parameters
# * Dependent variable: gamble
# * Get the RSS and R^2

# Independent variables
X = model.matrix( ~ sex + status + income + verbal + income:sex, mydata)
# Response or dependent variable
y = mydata$gamble

# Lets apply the least squares formula
# ( (X' X)^-1 * x  X'  * y
# ?solve -> is is used to solve the inverse of a matrix
beta = solve(t(X) %*% X) %*% t(X) %*% y

# Estimated output
yhat = X %*% beta

# Estimation Error
error = y - yhat

# Residual Sum of Squares
RSS = t(error) %*% error
sprintf("RSS = %f", as.double(RSS))

# R-squared
numerator = sum((yhat-y)^2)
denominator = sum((y - mean(yhat))^2)
R.squared = 1 - (numerator/denominator)
sprintf("R^2 = %f", as.double(R.squared))

# * Use the R function to compare your answers
# Compute the linear model using R functions
lmFitModel = lm(formula = gamble ~ sex + status + income + verbal + income:sex, data = mydata)
#print the model's summary to verify the R-squared calculation
summary(lmFitModel)
# Analysis of Variance Model to verify the RSS calculation
aov(lmFitModel)
#Let's remove sex, status, verbal and sex as they do not have a significant impact on the model
lmFitModel = lm(formula = gamble ~ income + income:sex, data = mydata)

# Which regression did you use?
# * ANCOVA, ANOVA, simple regression, logistic regression
# * Justify your answer
sprintf("ANCOVA, since the predictors are a mixture of quantitative and qualitative.")



