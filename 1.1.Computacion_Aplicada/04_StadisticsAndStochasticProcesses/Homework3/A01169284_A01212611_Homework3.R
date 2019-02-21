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
#*     - http://www.di.fc.ul.pt/~jpn/r/EM/GaussianMix.html
#*     - http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.409.3306&rep=rep1&type=pdf
#*
#* START DATE :
#*     18 Feb 2019
#*******************************************************************************

# Set worging directory to source file location
setwd("C:/Git/MNT_ITESM_courses/1.1.Computacion_Aplicada/04_StadisticsAndStochasticProcesses/Homework3");

# Install required libraries
# install.packages('e1071', dependencies=TRUE);
# install.packages('caret', dependencies=TRUE);
# install.packages('mixtools', dependencies=TRUE);
# install.packages('mclust', dependencies=TRUE);
# install.packages('BSDA', dependencies=TRUE);
# install.packages('VGAM', dependencies=TRUE);
# install.packages('psych', dependencies=TRUE);

# Load required libraries
library(caret);
library(datasets);
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

# From the problem we establish the Null hypothesis to be statistically the same
# in 2003 and 2004. The altern hypothesis is that the deaths in 2014 are
# significantly different than in 2003.

# Load the Dataset.csv
mydata = read.csv("Dataset.csv");

# complete.cases() returns a logical vector with the value TRUE for rows that
# are complete, and FALSE for rows that have some NA values
completeData = complete.cases(mydata);
# remove rows with incomplete data
mydata = mydata[completeData,];
head(mydata)
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
head(mydata)
# Only the counts related to "Drug use disorders"
# levels(mydata$Variable) #print the available values in column "Variable"
interestIndex = mydata$Variable == "Drug use disorders";
mydata        = mydata[interestIndex,];
head(mydata)
# Only the counts related to "Deaths per 100 000 population (standardised
# rates)"
# levels(mydata$Measure) #print the available values in column "Measure"
interestIndex = mydata$Measure == "Deaths per 100 000 population (standardised rates)";
mydata = mydata[interestIndex,];

# Let's take a look to the data ...
str(mydata);
summary(mydata);

# Only the counts related to "2003"
# levels(mydata$Year) #print the available values in column "Year"
interestIndex = mydata$Year == 2003;
my2003data    = mydata[interestIndex,];

# Only the counts related to "2014"
# levels(mydata$Year) #print the available values in column "Year"
interestIndex = mydata$Year == 2014;
my2014data    = mydata[interestIndex,];

# Get only relevant variables
# column | meaning
# -------+--------
#      3 | Country
#      5 | Value

columnsOfInterest = c(3,5);
# remove columns with repeated data
my2003data = my2003data[,columnsOfInterest];
my2014data = my2014data[,columnsOfInterest];

# When implementing a dependent t-test, both arguments shall have the same
# length, so let's remove the countries that dont appear in both datasets. (As
# we are dealing with paired samples)

Countries2Keep   =  my2003data$Country[my2003data$Country %in% my2014data$Country]
Countries2Remove =  my2003data$Country[!(my2003data$Country %in% my2014data$Country)]
# let's remove from my2003data: Canada France Italy Korea New Zealand
# Switzerland United Kingdom Slovenia Colombia 

summary(my2003data)
summary(my2014data)

# When implementing a dependent t-test, both arguments shall have the same length, so let's remove the countries that dont appear in both datasets. (As we are dealing with paired samples)
# Finding countries to remove from the 2003 set:
Countries2Keep   = my2003data$Country[my2003data$Country %in% my2014data$Country]
Countries2Keep
Countries2Remove = my2003data$Country[!(my2003data$Country %in% my2014data$Country)]
Countries2Remove

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

# Finding countries to remove from the 2014 set:
Countries2Keep   = my2014data$Country[my2014data$Country %in% my2003data$Country]
Countries2Remove = my2014data$Country[!(my2014data$Country %in% my2003data$Country)]
Countries2Remove

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
#   * the data set contains the same subjects (countries) measured on the same
#     variable twice (in 2003 and 2014).
#   * if we have the same country measured twice, then we calcualte a difference
#     score for each country,
#   * and then the mean of the difference scores.
#   * Notice that the t-value is the observed mean of the difference scores
#     relative to a standard error, which is the difference we expect by chance
#     due to sampling error.

# Performs one and two sample t-tests on vectors of data.
t.test(my2003data$Value,my2014data$Value,paired = TRUE);

# Where t = (Observed - Expected) / StdError
# Observed = Mean of difference values
# Expected = Zero "no effect", to assume the null hypotesis is true
# StdError = Standard error of the mean of the difference scores

# If the p-value is inferior or equal to 0.05, we can conclude that the
# difference between the two paired samples are significantly different.

Difference = my2014data$Value - my2003data$Value
plot(Difference,
     pch = 16,
     ylab="Difference (my2014data$Value - my2003data$Value)")
abline(0,0, col="#7DB0DD", lwd=2)
# A simple plot of differences between one sample and the other. Points below
# the blue line indicate observations where my2003data$Value is greater than
# my2014data$Value, that is where (my2014data$Value - my2003data$Value) is
# negative.

plot(my2003data$Value, my2014data$Value,
     pch = 16,
     xlab="my2003data$Value",
     ylab="my2014data$Value")
abline(0,1, col="#7DB0DD", lwd=2)
# Plot of paired samples from a paired t-test.  Circles below or to the right of
# the blue one-to-one line indicate observations with a higher value for
# my2003data$Value than for my2014data$Value.

# On the other hand, t-test can be biased by sample size ...
# When the sample size is very large, the standard error is very small - So
# even a small observed difference may be stadistically significant.

# Since our sample size is 24, the magnitude of the effect is not significant
# Therefore the t-test is not biased by the sample size.

x <- Difference
h <- hist(x, breaks=50, col="#7DB0DD", xlab="Difference", 
        main="Histogram of differences") 
xfit <- seq(min(x),max(x),length=24) 
yfit <- dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="#86B875", lwd=2)

# Histogram of differences of two populations from a paired t-test.  Distribution
# of differences should be approximately normal.  Bins with negative values
#indicate observations with a higher value for my2003data$Value than for my2014da
# a$Value.
# There is enough evidence to reject the Null Hypothesis, thus the deaths in 2014
# are significantly different than in 2003.


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
species_data <- mydata[,-5]
species_col <- c("#7DB0DD","#86B875","#E495A5")

# let's use "EM ALGORITHM FOR MIXURES OF MULTIVARIATE NORMALS" to cluster values
# into 3 classes

# Ignoring the known labels (species) of the Fisher Iris data, let us identify
# three clusters with the k-means method and compute the missclassification rate
set.seed(1234)

# labels are the original ones with this seed (avoid permutation)
r.km <- kmeans(species_data, centers=3)
mean(r.km$cluster!=as.numeric(mydata$Species))*100

# Let us know fit a mixture of three multidimensional Gaussian distributions.
# The model assumes the same variance covariance matrix for the three
# distributions (arbvar=FALSE). Initial centers are those given by the kmeans
# procedure.
library(mixtools)
c0 <- list(r.km$centers[1,], r.km$centers[2,], r.km$centers[3,])
mixmdl <- mvnormalmixEM(species_data, mu=c0, arbvar=FALSE)
summary(mixmdl) # lambda is the proportion of each cluster

# Fitting mixtures - get confusion matrix
pred <- apply(mixmdl$posterior, 1, function(row) which.max(row))
pred <- mapvalues(pred,
                  from = c("1", "2", "3"),
                  to   = c("setosa", "versicolor", "virginica"))
table(species_labels, pred)

# Let's plot

#CLUSTER PLOTS#
detach(package:mixtools)
plotClusters <- function(j1,j2){
  Species <- levels(mydata$Species)
  df_ell <- data.frame()
  for(g in (1:3)){
    M=mixmdl$sigma[c(j1,j2),c(j1,j2)]
    c=mixmdl$mu[[g]][c(j1,j2)]
    df_ell <- rbind(df_ell, cbind(as.data.frame(ellipse(
      M,
      centre=c,
      level=0.68)), group=Species[g]))
  }
  pl0 <- ggplot(data=mydata) +
    geom_point(aes_string(x=mydata[,j1],y=mydata[,j2], colour=species_labels))+ 
    theme(legend.position="bottom")+
    xlab(names(mydata)[j1])+
    ylab(names(mydata)[j2])+
    geom_path(data=df_ell, aes(x=x, y=y,color=group), size=1, linetype=1)
  return(pl0)
}
pl1 <- plotClusters(2,1) # Sepal.Width  vs Sepal.Length
pl2 <- plotClusters(3,1) # Petal.Length vs Sepal.Length
pl3 <- plotClusters(4,1) # Petal.Width  vs Sepal.Length
pl4 <- plotClusters(3,2) # Petal.Length vs Sepal.Width
pl5 <- plotClusters(4,2) # Petal.Width  vs Sepal.Width
pl6 <- plotClusters(4,3) # Petal.Width  vs Petal.Length
grid.arrange(pl1,pl2,pl3,pl4,pl5,pl6, layout_matrix=rbind(c(1,2,3),c(4,5,6)))

#DENSITY PLOTS#
plotDensities <- function(n){
  pl0 <- ggplot(mydata, aes(x=mydata[,n], color=Species, fill=Species))+
    geom_histogram(aes(y=..density..), alpha=0.5, position="identity")+
    geom_density(alpha=.2)+
    labs(title=names(mydata)[n],x="Value", y = "Density")
  return(pl0)
}
pl1 <- plotDensities(1)
pl2 <- plotDensities(2)
pl3 <- plotDensities(3)
pl4 <- plotDensities(4)
grid.arrange(pl1,pl2,pl3,pl4, layout_matrix=rbind(c(1,2),c(3,4)))

#BOX PLOTS#
plotBoxes <- function(n){
  pl0 <- ggplot(data=mydata, aes(x=Species, y=mydata[,n], color=Species))+
    geom_boxplot()+
    labs(title=names(mydata)[n],x="Species", y = "Value")+
    theme(legend.position="none")
  return(pl0)
}
pl1 <- plotBoxes(1)
pl2 <- plotBoxes(2)
pl3 <- plotBoxes(3)
pl4 <- plotBoxes(4)
grid.arrange(pl1,pl2,pl3,pl4, layout_matrix=rbind(c(1,2),c(3,4)))