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
#*     XX Feb 2019
#************************************************************************



#*** PART 1 *************************************************************

# Do a regression analysis to the “pima” dataset from the “faraway”
# library
library(faraway)
data(pima)

# * Analyze the database and select only the observations with no missing
#   data
# complete.cases() returns a logical vector with the value TRUE for rows
# that are complete, and FALSE for rows that have some NA values
completeData = complete.cases(pima)
# remove rows with incomplete data
mydata = pima[completeData, ]

# * Use the R function to fit a model.
#head(mydata) #returns: pregnant, glucose, diastolic, triceps, insulin,
#bmi, diabetes, age, test

# * Dependent variable : Test
# * How many correct predictions did the fitted model get? How many wrong?
# * Confusion matrix
# Which regression did you use?
# * ANCOVA, ANOVA, simple regression, logistic regression
# * Justify your answer



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


