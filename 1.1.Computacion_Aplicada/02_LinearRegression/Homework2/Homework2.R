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

# Do a regression analysis to the “pima” dataset from the “faraway”
# library
library(faraway)
library(caret)
data(pima)

# * Analyze the database and select only the observations with no missing
#   data
# complete.cases() returns a logical vector with the value TRUE for rows
# that are complete, and FALSE for rows that have some NA values
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

#------------------------------------------------------------------------

#LINEAR ANALYSIS
# lmFitModel <- lm(test ~ pregnant + glucose + diastolic + triceps + insulin + bmi + diabetes + age, data = mydata)
# summary(lmFitModel)
# lmFitModel
#Let's remove triceps, insulin, and age as they do not have a significant impact on the model
lmFitModel <- lm(test ~ pregnant + glucose + diastolic + bmi + diabetes, data = mydata)
summary(lmFitModel)
lmFitModel

#------------------------------------------------------------------------

#POLYNOMIAL ANALYSIS
gp6Model = lm(
  test ~
      (pregnant + glucose + diastolic + bmi + diabetes)    +
    I((pregnant + glucose + diastolic + bmi + diabetes)^2) +
    I((pregnant + glucose + diastolic + bmi + diabetes)^3) +
    I((pregnant + glucose + diastolic + bmi + diabetes)^4) +
    I((pregnant + glucose + diastolic + bmi + diabetes)^5) +
    I((pregnant + glucose + diastolic + bmi + diabetes)^6),
  data = mydata)
summary(gp6Model)
gp6Model

#------------------------------------------------------------------------

#ANOVA ANALYSIS
# aovFitModel <- aov(test ~ pregnant + glucose + diastolic + triceps + insulin + bmi + diabetes + age, data = mydata)
# summary(aovFitModel)
# aovFitModel
#Let's remove diastolic, triceps, and insulin as they do not have a significant impact on the model
aovFitModel <- aov(test ~ pregnant + glucose + bmi + diabetes + age, data = mydata)
summary(aovFitModel)
aovFitModel

#------------------------------------------------------------------------

# * How many correct predictions did the fitted model get? How many wrong?

#get test predictions using the LINEAR ANALYSIS
#We assume that: all values less than 0.5 correspond to the first category,
#0.5 or bigger correspond to the second category
lm_predictions = ifelse(predict(lmFitModel, mydata) > 0.5, "1", "0" )
#compare the actual data agains the predictions
correct_lm_predictions = mydata$test == lm_predictions
#count how many of the predictions are right vs wrong
#FALSE:WRONG TRUE:RIGHT
summary(correct_lm_predictions)
# * Confusion matrix
confusionMatrix(as.factor(lm_predictions), as.factor(mydata$test), positive = NULL, dnn = c("Prediction", "Reference"))

#get test predictions using the POLYNOMIAL ANALYSIS
#We assume that: all values less than 0.5 correspond to the first category,
#0.5 or bigger correspond to the second category
gp6_predictions = ifelse(predict(gp6Model, mydata) > 0.5, "1", "0" )
#compare the actual data agains the predictions
correct_gp6_predictions = mydata$test == gp6_predictions
#count how many of the predictions are right vs wrong
#FALSE:WRONG TRUE:RIGHT
summary(correct_gp6_predictions)
# * Confusion matrix
confusionMatrix(as.factor(gp6_predictions), as.factor(mydata$test), positive = NULL, dnn = c("Prediction", "Reference"))

#get test predictions using the ANOVA ANALYSIS
#We assume that: all values less than 0.5 correspond to the first category,
#0.5 or bigger correspond to the second category
aov_predictions = ifelse(predict(aovFitModel, mydata) > 0.5, "1", "0" )
#compare the actual data agains the predictions
correct_aov_predictions = mydata$test == aov_predictions
#count how many of the predictions are right vs wrong
#FALSE:WRONG TRUE:RIGHT
summary(correct_aov_predictions)
# * Confusion matrix
confusionMatrix(as.factor(aov_predictions), as.factor(mydata$test), positive = NULL, dnn = c("Prediction", "Reference"))

# Which regression did you use?
# * ANCOVA, ANOVA, simple regression, logistic regression
# * Justify your answer

#TODO: explain ...

#*** PART 2 *************************************************************

# Do a regression analysis to the “teengamb” dataset from the “faraway”
# library
data(teengamb)
# complete.cases() returns a logical vector with the value TRUE for rows
# that are complete, and FALSE for rows that have some NA values
completeData = complete.cases(teengamb)
# remove rows with incomplete data
mydata = teengamb[completeData, ]

# * Use the normal equations to fit the 𝛽 parameters
# * Dependent variable: gamble
# * Get the RSS and 𝑅2
# * Use the R function to compare your answers
# Which regression did you use?
# * ANCOVA, ANOVA, simple regression, logistic regression
# * Justify your answer


