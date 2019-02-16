#*******************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     Homework3.R
#*
#* DESCRIPTION :
#*     Computación Aplicada (Ene 19 Gpo 1)
#*     Homework 3
#*
#* NOTES :
#*     - 
#*
#* START DATE :
#*     12 Feb 2019
#*******************************************************************************

# Set worging directory to source file location
setwd("~/MNT_ITESM_courses/1.1.Computacion_Aplicada/04_StadisticsAndStochasticProcesses/Homework3");

# Install required libraries
#install.packages('e1071', dependencies=TRUE);
#install.packages('caret', dependencies=TRUE);
#install.packages('mixtools', dependencies=TRUE);
#install.packages('mclust', dependencies=TRUE);
#install.packages('BSDA', dependencies=TRUE);

# Load required libraries
library(caret);
library(datasets);
library(mixtools);
library(mclust);
library(ggplot2);
library(dplyr);

#*******************************************************************************
#* PART 1 **********************************************************************
# With the Dataset.csv, filtered by "Drug use disorders" and "Deaths per
# 100,000 population (standardized rates)" apply a statistical test to see if
# the deaths in 2014 are significantly different than in 2003. (50%)
# * Justify the answer and the use of the statistical test

# Load the Dataset.csv
mydata = read.csv("Dataset.csv");

# complete.cases() returns a logical vector with the value TRUE for rows that
# are complete, and FALSE for rows that have some NA values
completeData = complete.cases(mydata);
# remove rows with incomplete data
mydata = mydata[completeData,];

# Get only relevant variables (Not repeated data accross columns)
# column | meaning
# -------+--------
#      2 | Variable
#      4 | Measure
#      6 | Country
#      8 | Year
#      9 | Value
columnsOfInterest = c(2,4,6,8,9);
# remove columns with repeated data
mydata = mydata[,columnsOfInterest];

# Only the counts related to "Drug use disorders"
# levels(mydata$Variable) #print the available values in column "Variable"
interestIndex = mydata$Variable == "Drug use disorders";
mydata        = mydata[interestIndex,];

# Only the counts related to "Deaths per 100 000 population (standardised
# rates)"
# levels(mydata$Measure) #print the available values in column "Measure"
interestIndex = mydata$Measure == "Deaths per 100 000 population (standardised rates)";
mydata        = mydata[interestIndex,];

# Only the counts related to "2003"
# levels(mydata$Year) #print the available values in column "Year"
interestIndex = mydata$Year == 2003;
my2003data    = mydata[interestIndex,];

# Only the counts related to "2014"
# levels(mydata$Year) #print the available values in column "Year"
interestIndex = mydata$Year == 2014;
my2014data    = mydata[interestIndex,];

# Let's take a look to the data ...
str(mydata);
summary(mydata);

# Get only relevant variables
# column | meaning
# -------+--------
#      4 | Year
#      5 | Value
columnsOfInterest = c(4,5);
# remove columns with repeated data
my2003data = my2003data[,columnsOfInterest];
my2014data = my2014data[,columnsOfInterest];

# Let's take a look to the data ...
str(my2003data);
summary(my2003data);

str(my2014data);
summary(my2014data);

## Statistical tests
# Performs one and two sample t-tests on vectors of data.
t.test(my2003data$Value,my2014data$Value);
cat(sprintf("
mean of my2003data$Value: 0.4727273 deaths
mean of my2014data$Value: 0.2769231 deaths
    per 100,000 population (standardised rates) due to Drug use disorders.
\n"));

# Calculate the percentage decrease
decrease = 0.4727273 - 0.2769231;
decrease = decrease*100 / 0.4727273;

# Conclusion
cat(sprintf("
The deaths in 2014 are NOT significantly different than in 2003 as they only
decreased by 41.42 percent - (less than 50 percent).

The statistical test implemented was \"T-test (with 2 Samples)\" as it
contrasts the mean of two populations.
\n"));

# Performs one- and two-sample Wilcoxon tests on vectors of data; the latter is
# also known as 'Mann-Whitney' test.
#wilcox.test(my2003data$Value,my2014data$Value)
# [WARNING] message:
#   In wilcox.test.default(my2003data$Value, my2014data$Value) :
#   cannot compute exact p-value with ties

# Test for association between paired samples, using one of Pearson's product
# moment correlation coefficient, Kendall's tau or Spearman's rho.
#cor.test(my2003data$Value,my2014data$Value)
# [ERROR] in cor.test.default(my2003data$Value, my2014data$Value) : 
#   'x' and 'y' must have the same length

#*******************************************************************************
#* PART 2 **********************************************************************
# Dataset: iris
# * It gives the measurements in centimeters of the variables sepal length and
#   width and petal length and width, respectively, for 50 flowers from each of
#   3 species of iris.
# * Use a GMM algorithm to cluster values into 3 classes
# * Independent variables: Sepal length and width, petal length and width
# * Return Confusion matrix of predicted class versus real class (Species)
# * Plot the answer (+10)
# * Independent variables or transformation of independent variables colored by
#   the predicted class

# Load the iris dataset
data(iris)

# complete.cases() returns a logical vector with the value TRUE for rows that
# are complete, and FALSE for rows that have some NA values
completeData = complete.cases(iris);
# remove rows with incomplete data
mydata = iris[completeData,];

# Let's take a look to the data ...
str(mydata);
summary(mydata);
pairs(mydata);
#cor(mydata);

###########################

mb = Mclust(mydata[,-5])

#or specify number of clusters
mb3 = Mclust(mydata[,-5], 3)

# optimal selected model
mb3$modelName

# optimal number of cluster
mb3$G

# probality for an observation to be in a given cluster
head(mb3$z)

# get probabilities, means, variances
summary(mb3, parameters = TRUE)

# Compare amount of the data within each cluster
table(iris$Species, mb$classification)
# vs
table(iris$Species, mb3$classification)

# After the data is fit into the model, we plot the model based on clustering results.
plot(mb3, "classification")
# The types of what argument are: "BIC", "classification", "uncertainty", "density". By default, all the above graphs are produced.

# Let's plot estimated density. Contour plot
plot(mb3, "density")
plot(mb3, "uncertainty")

# You can also use the summary() function to obtain the most likely model and the most possible number of clusters.
summary(mb3)

###########################

# #There is data from 3 Species
# #Let's convert the Species data into a factor
# # 1:setosa
# # 2:versicolor
# # 3:virginica
# mydata$Species = factor(as.numeric(mydata$Species))
# 
# measurements =
#   mydata$Sepal.Length +
#   mydata$Sepal.Width +
#   mydata$Petal.Length +
#   mydata$Petal.Width;
# hist(measurements)
# 
# # Gaussian Mixture Model (GMM)
# # normalmixEM() K indica la cantidad de estados ocultos que queremos usar
# # EM indica el uso del algoritmo de esperanza maximización
# mixmdl <- normalmixEM(measurements,k=3);
# summary(mixmdl); # lambda significa la proporción que hay en cada grupo
# 
# # mixtools sobreescribe la función plot
# # Con el valor 2 indicamos que queremos imprimir las curvas de densidad
# plot(mixmdl,whichplots=2) 
# lines(density(measurements), lty=2, lwd=2)
# 
# # Análisis bidimensional
# colBind = as.data.frame(cbind(measurements,mydata$Species))
# names(colBind) <- c("measurements","Species")
# 
# plot(colBind)
# mixmv <- mvnormalmixEM(colBind,k=3)
# summary(mixmv)
# plot(mixmv,which=2)

