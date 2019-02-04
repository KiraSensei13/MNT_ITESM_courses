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