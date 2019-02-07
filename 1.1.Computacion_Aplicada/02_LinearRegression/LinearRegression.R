# Install package/library 
# install.packages("faraway")
library(faraway)
data(pima)

# Exploration
head(pima)
summary(pima)

# Something weird out there?
#bmi, glucosa, triceps, insulin, diastolic != 0
# Test -> categoric data
plot(diabetes ~ as.factor(test), pima)

############################################################ Lecture03 follows below

# R = R

# Load library
library(faraway)
data(gala)

#Exploration
head(gala)
dim(gala)

# How do we set our matrix in R?
# Independent variables
X = model.matrix( ~ Endemics + Area + Elevation + Nearest + Scruz + Adjacent , data=gala )
# Response or dependent variable
y = gala$Species

# Lets apply the least squares formula
# ( (X' X)^-1 * x  X'  * y
# ?solve -> is is used to solve the inverse of a matrix
beta = solve(t(X) %*% X) %*% t(X) %*% y
beta

# Estimated output
yhat = X %*% beta

# Estimation Error
error = y - yhat

# Residual Sum of Squares
RSS = t(error) %*% error

# R-squared
numerator = sum((yhat-y)^2)
denominator = sum((y - mean(yhat))^2)
R.squared = 1 - (numerator/denominator)

# Compute the linear model using R functions
mdl = lm( Species ~ Endemics + Area + Elevation + Nearest + Scruz + Adjacent , data=gala)
summary(mdl)

## Lets ask for the rank of a matrix
library("Matrix")

rankMatrix(gala)
rankMatrix(gala)[1]

# Lets add a new variable
diff = gala$Nearest - gala$Scruz
gala$diff = diff

mdl = lm( Species ~ Endemics + Area + Elevation + Nearest + Scruz + Adjacent + diff, data=gala)
summary(mdl)

# Did the rank changed?
rankMatrix(gala)[1]

# Lets add a white noise
gala$diff = diff + rnorm( length(diff), sd = 0.0001 )

mdl = lm( Species ~ Endemics + Area + Elevation + Nearest + Scruz + Adjacent + diff, data=gala)
summary(mdl)

# Did the rank changed?
rankMatrix(gala)[1]

# Second part
## Polynomial regression
library(faraway)
data("corrosion")
head(corrosion)

# It is necessary to use the I() function to use the polynomial regression

gp2 = lm( loss ~ Fe + I(Fe^2), corrosion )
gp3 = lm( loss ~ Fe + I(Fe^2)+ I(Fe^3) , corrosion )
gp4 = lm( loss ~ Fe + I(Fe^2)+ I(Fe^3)+ I(Fe^4) , corrosion )
gp5 = lm( loss ~ Fe + I(Fe^2)+ I(Fe^3)+ I(Fe^4)+ I(Fe^5) , corrosion )
gp6 = lm( loss ~ Fe + I(Fe^2)+ I(Fe^3)+ I(Fe^4)+ I(Fe^5)+ I(Fe^6) , corrosion )


# How can we plot our answers??
# Lets use ?predict

# predict(  lm estimated with R, data.frame with the new values to predict   (X matrix) )
grid = seq(0,2,len = 50)

# seq is a function used to generate a sequence of values
#seq ( initial value, last value, step bewtween generated values)

plot(loss ~ Fe, data = corrosion, ylim=c(60, 140))
lines(grid, predict(gp2, data.frame(Fe=grid)) , col = "red")

plot(loss ~ Fe, data = corrosion, ylim=c(60, 140))
lines(grid, predict(gp3, data.frame(Fe=grid)) , col = "blue")

plot(loss ~ Fe, data = corrosion, ylim=c(60, 140))
lines(grid, predict(gp4, data.frame(Fe=grid)) , col = "green")

plot(loss ~ Fe, data = corrosion, ylim=c(60, 140))
lines(grid, predict(gp5, data.frame(Fe=grid)) , col = "magenta")

plot(loss ~ Fe, data = corrosion, ylim=c(60, 140))
lines(grid, predict(gp6, data.frame(Fe=grid)) , col = "chocolate")

# El mejor parece ser un polinomio de grado 6
plot(loss ~ Fe, data = corrosion, ylim=c(60, 130))
points(corrosion$Fe, fitted(gp6), pch = 18, col = "chocolate")
lines(grid, predict(gp6, data.frame(Fe=grid)) , col = "chocolate")

#####################################################################

# R = R + N

#R:numerical  N:categorical

## ANCOVA
data(sexab)
head(sexab)

plot(ptsd~csa, sexab)
plot(ptsd~cpa, pch=as.character(csa), sexab)

# If the lm functions detects a text , it automatically changes the response to a categorical value
g = lm(ptsd ~ cpa+csa+cpa:csa, sexab)
summary(g)

g = lm(ptsd ~ cpa+csa, sexab)
summary(g)

# The base slope is 0.5506 (cpa)
# The second slope is -6.273 (csaNotAbused)
# Lets plot them!
plot(ptsd~cpa, pch=as.character(csa), sexab)
abline(10.2480, 0.5506, col = "red")
abline(10.2480-6.2728, 0.5506, col = "magenta")

#####################################################################

# R = N

## ANOVA
data("coagulation")
head(coagulation)

plot(coag ~ diet, coagulation, ylab="Coagulation time")

# Assumption mu = 0
g2 = lm(coag ~ diet - 1, coagulation)
summary(g2)
model.matrix(g2)

# Assumption alpha 1 = 0, where alpha 1 is the diet1 beta value
g = lm(coag ~ diet, coagulation)
summary(g)
model.matrix(g)

#####################################################################

# N = R

# If I set alpha1 to zero, I cannot use the R-squared value
## Logistic Regression
#  R function to use: glm
lrm = glm(csa ~ cpa + ptsd, data=sexab, family=binomial)

# Prediction to be ploted. As we use a logistic model we have to specify the output type as response
prediction = predict(lrm, data.frame(cpa=sexab$cpa, ptsd = sexab$ptsd), type = "response") 
hist(prediction)

# What happened? I need 0 or 1, and I have a lot more!!!
# We assume that: all values less than 0.5 correspond to the first category, whilst 0.5 or bigger correspond to the second condition
# How to do it??
newPrediction = ifelse(prediction > 0.5, "NotAbused", "Abused" )

# Confussion matrix to see the model performance
variableTable = table(newPrediction , as.character(sexab$csa))

# get the percentage of prediction ...
sum(diag(variableTable))/sum(variableTable)




table(newPrediction , as.character(sexab$csa))
