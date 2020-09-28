## Load Library
library(rpart)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(knitr)
library(lattice)
library(dplyr)



## set wd
setwd("~/SMU MS DS/MSDS 7333/wk5")


## Load Rda 
load(file = "~/SMU MS DS/MSDS 7333/wk5/data.rda")
#head(emailDFrp)
## extra step not required, this is to get a dataset as CSV
#write.csv(emailDFrp,file="data.csv")


## Data Reading
## Section 3.3

## list directories
colnames(emailDFrp)

## Finding the total length of the spam dataset
lengths(emailDFrp)


## splitting Test and Train
# total numbers of email in list, emailDFrp[1] is the isSpam column
numEmail<- as.numeric(lengths(emailDFrp[1]))
# total number of Spam
numSpam <- as.numeric(sum(emailDFrp[1]!="T"))
# total number of Ham
numHam <- numEmail - numSpam

# set random seed for repeat selection process
set.seed(123456)

# Spam and Ham test indices
testSpamIdx <- sample(numSpam, size = floor(numSpam/3))
testHamIdx <- sample(numHam, size = floor(numHam/3))

#msgWordsList function?


