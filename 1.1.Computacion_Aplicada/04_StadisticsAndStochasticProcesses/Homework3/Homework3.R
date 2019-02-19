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
#*     - http://sia.webpopix.org/mixtureModels.html#the-iris-data
#*     - https://cran.r-project.org/web/packages/mclust/vignettes/mclust.html
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
#install.packages('VGAM', dependencies=TRUE);
#install.packages('psych', dependencies=TRUE);

# Load required libraries
library(caret);
library(datasets);
library(mixtools);
library(mclust);
library(ggplot2);
library(dplyr);
library(VGAM);
library(ellipse);
library(gridExtra);
library(psych);

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
#      3 | Country
#      5 | Value
columnsOfInterest = c(3,5);
# remove columns with repeated data
my2003data = my2003data[,columnsOfInterest];
my2014data = my2014data[,columnsOfInterest];

# When implementing a dependent t-test, both arguments shall have the same length, so let's remove the countries that dont appear in both datasets. (As we are dealing with paired samples)

Countries2Keep   = my2003data$Country[my2003data$Country %in% my2014data$Country]
Countries2Remove = my2003data$Country[!(my2003data$Country %in% my2014data$Country)]
# let's remove from my2003data: Canada France Italy Korea New Zealand Switzerland United Kingdom Slovenia Colombia 
interestIndex = my2003data$Country != "Canada" &
                my2003data$Country != "France" &
                my2003data$Country != "Italy" &
                my2003data$Country != "Korea" &
                my2003data$Country != "New Zealand" &
                my2003data$Country != "Switzerland" &
                my2003data$Country != "United Kingdom" &
                my2003data$Country != "Slovenia" &
                my2003data$Country != "Colombia";
my2003data    = my2003data[interestIndex,];

Countries2Keep   = my2014data$Country[my2014data$Country %in% my2003data$Country]
Countries2Remove = my2014data$Country[!(my2014data$Country %in% my2003data$Country)]
# let's remove from my2014data: Slovak Republic Estonia
interestIndex = my2014data$Country != "Slovak Republic" &
                my2014data$Country != "Estonia";
my2014data    = my2014data[interestIndex,];

# Let's take a look to the data ...
str(my2003data);
summary(my2003data);

str(my2014data);
summary(my2014data);

## Statistical tests

# Dependent t-test is to be implemented ...
# AKA. paired samples t-test - the reasons are the following:
#   * the data set contains the same subjects (countries) measured on the same variable twice (in 2003 and 2014).
#   * if we have the same country measured twice, then we calcualte a difference score for each country,
#   * and then the mean of the difference scores.
#   * Notice that the t-value is the observed mean of the difference scores relative to a standard error, which is the difference we expect by chance due to sampling error.

# Performs one and two sample t-tests on vectors of data.
t.test(my2003data$Value,my2014data$Value,paired = TRUE);

# Where t = (Observed - Expected) / StdError
# Observed = Mean of difference values
# Expected = Zero "no effect", to assume the null hypotesis is true
# StdError = Standard error of the mean of the difference scores

# If the p-value is inferior or equal to 0.05, we can conclude that the difference between the two paired samples are significantly different.

Difference = my2014data$Value - my2003data$Value
plot(Difference,
     pch = 16,
     ylab="Difference (my2014data$Value - my2003data$Value)")
abline(0,0, col="#7DB0DD", lwd=2)
# A simple plot of differences between one sample and the other.  Points below the blue line indicate observations where my2003data$Value is greater than my2014data$Value, that is where (my2014data$Value - my2003data$Value) is negative.

plot(my2003data$Value, my2014data$Value,
     pch = 16,
     xlab="my2003data$Value",
     ylab="my2014data$Value")
abline(0,1, col="#7DB0DD", lwd=2)
# Plot of paired samples from a paired t-test.  Circles below or to the right of the blue one-to-one line indicate observations with a higher value for my2003data$Value than for my2014data$Value.

# On the other hand, t-test can be biased by sample size ...
#   when the sample size is very large, the standard error is very small - So even a small observed difference may be stadistically significant

# Since our sample size is 24, the magnitude of the effect is not significant - and therefore the t-test is not biased by the sample size

x <- Difference
h<-hist(x, breaks=50, col="#7DB0DD", xlab="Difference", 
        main="Histogram of differences") 
xfit<-seq(min(x),max(x),length=40) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="#86B875", lwd=2)
# Histogram of differences of two populations from a paired t-test.  Distribution of differences should be approximately normal.  Bins with negative values indicate observations with a higher value for my2003data$Value than for my2014data$Value.


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

species_labels <- mydata[,5]
species_col <- c("#7DB0DD","#86B875","#E495A5")

pairs(mydata[,-5], col = species_col[species_labels],
      lower.panel = NULL, cex.labels=2, pch=19, cex = 1)
# "setosa"BLUE "versicolor"GREEN "virginica"RED
par(xpd = TRUE)
legend(
  "bottomleft",
  inset=0.02,
  legend = as.character(levels(species_labels)),
  fill = species_col)

# create scatterplots, histograms & correlation coefficients
pairs.panels(mydata[,-5],
             gap=0,
             bg=species_col[mydata$Species],
             pch=21)

summary(mydata)

# #Data normalization and distance matrix
# # variables that have larger values, they contribute more, and the entire clustering will be dominated by them. So to avoid that, what we do is we normalize all the variables so that the average is 0 and standard deviation is 1, so let's do that
# z = mydata[,-5]
# m = apply(z,2,mean)
# s = apply(z,2,sd)
# 
# z = scale(z,center=m, scale=s)
# mydata = cbind(data.frame(z), mydata[5]) # now all values are between -3 and 3
# 
# pairs.panels(mydata[,-5],
#              gap=0,
#              bg=species_col[mydata$Species],
#              pch=21)
# 
# summary(mydata)

#Gaussian clustering for iris data#

# ...












