### HW 2

## Part 1
library(faraway)
data(pima)

# Exploration
summary(pima)

# Apply filter by filter
index = pima$glucose > 0
pima = pima[ index,  ]

index = pima$diastolic > 0
pima = pima[ index,  ]

index = pima$insulin > 0
pima = pima[ index,  ]

index = pima$triceps > 0
pima = pima[ index,  ]

index = pima$bmi > 0
pima = pima[ index,  ]

filteredPima = pima

# # All the filters at the same time
# index <- with(pima ,
#                   glucose > 0 &
#                   diastolic > 0 &
#                   triceps > 0 &
#                   insulin > 0 &
#                   bmi > 0
#                 )
# filteredPima <- pima[ index , ]

# New size
dim(filteredPima) # 392 x 9
summary(filteredPima)


# Set numeric values to categoric value
filteredPima$test <- as.factor(filteredPima$test)
summary(filteredPima)

# Which linear model am I going to use?
# Logistic regression
mdl <- glm( test ~ ., data = filteredPima, family = "binomial")
summary(mdl)

mdl <- glm( test ~ glucose + bmi + diabetes + age, data = filteredPima, family = "binomial")
summary(mdl)


# Lets see the performance of our beta coefficients
predicted <- predict(mdl, filteredPima[,-9], type = "response")
hist(predicted)
predicted <- ifelse(predicted < 0.5, "0", "1")

predicted.table = table("Real" = filteredPima[,9], "Predicted" = predicted)
print(predicted.table)
# Prediction Power
print(sum(diag(predicted.table))/sum(predicted.table))

## Part 2
data("teengamb")
head(teengamb)
summary(teengamb)

# ANCOVA regression
teengamb$sex <- as.factor(teengamb$sex)
summary(teengamb)

# Lets create our inependent variable matrix and dependent variable vector
X <- model.matrix( ~  sex + status + income + verbal, data=teengamb )
y <- teengamb$gamble

# First formula
beta <- solve( t(X) %*% X) %*% t(X) %*% y 

# Estimates
yhat <- X %*% beta
error <- y - yhat

# Residual Sum of Squares
RSS <- t(error) %*% error

# R-squared
numerator = sum( ( yhat - y)^2 )
denominator = sum( ( y - mean(y))^2 )

R.squared = 1 - numerator/denominator


# Using R function
mdl <- lm( gamble ~ . , data=teengamb)
summary(mdl)

# Lets compare
print(beta)
print(RSS)
print(R.squared)
