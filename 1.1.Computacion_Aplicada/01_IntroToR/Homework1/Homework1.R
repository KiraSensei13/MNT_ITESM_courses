#************************************************************************
#* AUTHOR(S) :
#*     Bruno González Soria          (A01169284)
#*     Antonio Osamu Katagiri Tanaka (A01212611)
#*
#* FILENAME :
#*     Homework1.R
#*
#* DESCRIPTION :
#*     Computación aplicada (Ene 19 Gpo 1)
#*     Homework 1
#*
#* NOTES :
#*     - Upload the .R file with all your code in the Assignment tab in
#*       Bb Assignment - R Homework
#*     - 1 per pair. Both names in a comment at the beginning of the .R
#*       file
#*     - DO NOT upload the csv file
#*
#* START DATE :
#*     25 Jan 2019
#************************************************************************



#*** PART 1 *************************************************************

# Load the dataset.csv
dataset = read.csv("Dataset.csv")

# The new dataset must have:
# * Only relevant variables (Not repeated data accross columns)
# column | meaning
# -------+--------
#      2 | Variable
#      4 | Measure
#      6 | Country
#      8 | Year
#      9 | Value
columnsOfInterest = c(2,4,6,8,9)
mydataset         = dataset[, columnsOfInterest]

# * Only the counts related to "All causes of death"
# levels(mydataset$Variable) #print the available values in column "Variable"
interestIndex = mydataset$Variable == "All causes of death"
mydataset     = mydataset[interestIndex, ]

# * And the measure must be "Number of total deaths"
# levels(mydataset$Measure) #print the available values in column "Measure"
interestIndex = mydataset$Measure == "Number of total deaths"
mydataset     = mydataset[interestIndex, ]

# * All years and all countries must be included
# * Save the new database in another csv file called "filteredDatabase.csv"
# TIP: Use matrix and submatrix notation to get the desired values and the write.csv function
write.csv(mydataset, file = "filteredDatabase.csv", row.names = FALSE, na = "")
#     remove row names, omit NAs



#*** PART 2 *************************************************************

# Select another cause of death and measure to get a second dataset.
# Filters to select the number of male deaths due to VIH-AIDS
mydataset     = dataset[, columnsOfInterest]
interestIndex = mydataset$Variable == "HIV-AIDS"
mydataset     = mydataset[interestIndex, ]
interestIndex = mydataset$Measure == "Number of male deaths"
mydataset     = mydataset[interestIndex, ]

# With this new data set:
# * Plot the boxplot of these deaths in the first year reported and compare it with the last year
#   (both boxplots in the same plot)
FirstYear     = min(mydataset$Year, na.rm = TRUE) # get the 1st year, remove NAs
LastYear      = max(mydataset$Year, na.rm = TRUE) # get the last year, remove NAs
# Apply filter
interestIndex = mydataset$Year == FirstYear | mydataset$Year == LastYear
mydataset     = mydataset[interestIndex, ]
#boxplot w/no outliners
boxplot(
  Value ~ Year,
  data = mydataset,
  xlab = "Year",
  ylab = "number of male deaths due to VIH-AIDS",
  outline = FALSE
)


