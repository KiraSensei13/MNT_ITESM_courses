### Data Analysis in R ###

#read data
dataset=read.csv("Dataset.csv")

#filter the columns that have repeated data
columnsOfInterest=c(2,4,6,8,9)
filteredDataset=dataset[,columnsOfInterest]

#1st filter to fetch rows where Variable="All causes of death"
#levels(filteredDataset$Variable) #print the available values in column "Variable"
interestIndex = filteredDataset$Variable == "All causes of death"
filteredDataset = filteredDataset[interestIndex,]

#2nd filter to fetch rows where Measure="Number of total deaths"
#levels(filteredDataset$Measure) #print the available values in column "Measure"
interestIndex = filteredDataset$Measure == "Number of total deaths"
filteredDataset = filteredDataset[interestIndex,]

#create a Box-plot with all the years
boxplot(Value ~ Year, data = filteredDataset, outline=F) #boxplot w/no outliners

#3rd filter to fetch data of 2010 & 2012 only
interestIndex = filteredDataset$Year == 2010 | filteredDataset$Year == 2012
filteredDataset = filteredDataset[interestIndex,]

#create a Box-plot for 2010 and 2012 only
#Value & Year are columns in Dataset.csv
boxplot(Value ~ Year, data = filteredDataset) #standard boxplot
boxplot(Value ~ Year, data = filteredDataset, outline=F) #boxplot w/no outliners
