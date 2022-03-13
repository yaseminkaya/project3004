#Load data

#datafolder <- "C:/Users/yyase/Downloads/Core Project Data/"
datafolder <- "C:/Users/Punkt/Downloads/Core Project Data/"

deaths <- read.csv(paste0(datafolder, "BBS3004_deaths.csv"), header = TRUE)
demo <- read.csv(paste0(datafolder, "BBS3004_demographics.csv"), header = TRUE)
hosp <- read.csv(paste0(datafolder, "BBS3004_hospitalisations.csv"), header = TRUE)
lab <- read.csv(paste0(datafolder, "BBS3004_labvalues.csv"), header = TRUE)
visits <- read.csv(paste0(datafolder, "BBS3004_visits.csv"), header = TRUE)


#Invert the table of labvalues #pivot_wider
#.libPaths("C:/Users/Punkt/Downloads/R/RStudio")

library(magrittr)
library(tidyr)

lab_inv <- lab %>%
  pivot_wider(names_from = biomarker, values_from = value)


#Create a common table for visits, labs and demo #reason for hospitalization

merged_table <- merge(x = visits, y = merge(x = merged_table, y = demo, by = "id"))
merged_table = merged_table[,!grepl(".x$",names(merged_table))]
merged_table = merged_table[,!grepl(".y$",names(merged_table))]


#initialize to 0, set back to 1
#hospitalization within 60

for (a in 1:nrow(merged_table)) {
  for (b in 1:nrow(hosp)) {
    if (merged_table$id[a] == hosp$id[b]) {
      if (difftime((as.Date(hosp$date[b])),
                   (as.Date(merged_table$date[a])),
                   units = "days") <= 60 &
          difftime((as.Date(hosp$date[b])),
                   (as.Date(merged_table$date[a])), 
                   units = "days") > 0) {
        merged_table$label_hosp[a] <- '1'
        break
      }
      else {
        merged_table$label_hosp[a] <- '0'
      }
    }
  }
}


#Death within 60 days

for (a in 1:nrow(merged_table)) {
  for (c in 1:nrow(deaths)) {
    if (merged_table$id[a] == deaths$id[c]) {
      if (is.na(deaths$date_of_death[c])) {
        merged_table$label_death[a] <- '0'
      }
      else {
        if (difftime((as.Date(deaths$date_of_death[c])),
                     (as.Date(merged_table$date[a])),
                     units = "days") <= 60 &
            difftime((as.Date(deaths$date_of_death[c])),
                     (as.Date(merged_table$date[a])), 
                     units = "days") > 0) {
          merged_table$label_death[a] <- '1'
          break
        }
        else {
          merged_table$label_death[a] <- '0'
        }
      }
    }
  }
}


#Death or hospitalization in 60 days

for (a in 1:nrow(merged_table)) {
  merged_table$label[a] <- 
    if (merged_table$label_hosp[a] == 1 || merged_table$label_death[a] == 1){
      1
    } else {0}}
merged_table = merged_table[,!grepl("hosp$",names(merged_table))]
merged_table = merged_table[,!grepl("death$",names(merged_table))]


#Outliers IQR for numerical predictors
#No outliers for Age and Biomarkers, except NT-BNP (93 observations)

outliers <- function(x) {
  
  Q1 <- quantile(x, .25)
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  
  outliers <- subset(merged_table, x < (Q1 - 1.5*IQR) | x > (Q3 + 1.5*IQR))
  
  return(outliers)
}

A <- outliers(merged_table$`NT-BNP`)
hist(log(merged_table$`NT-BNP`))
#outliers in NT.BNP --> not normally distributed


#Update age of the patient

merged_table$data_diff<-difftime(as.Date(merged_table$date), as.Date(merged_table$date_of_diagnosis), units = "days")
merged_table$age <- merged_table$age + (merged_table$data_diff/365)
merged_table = merged_table[,!grepl("^data_diff",names(merged_table))]


#Imputation of missing values - MICE

#Install necessary packages
install.packages("mice")
install.packages("VIM")
install.packages("Rcpp")
library("mice")
library("VIM")
library("Rcpp")


#Visualize what is missing

pMiss <- function(x){sum(is.na(x))/length(x)*100}
apply(merged_table,2,pMiss)
apply(merged_table,1,pMiss)
md.pattern(merged_table)
aggr_plot <- aggr(merged_table, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(data), cex.axis=.7, gap=3, 
                  ylab=c("Histogram of missing data","Pattern"))


#Necessary columns as factor, can also be already done in earlier part of code

merged_table$orthopnea <- as.factor(merged_table$orthopnea)
merged_table$oedema <- as.factor(merged_table$oedema)
merged_table$cough <- as.factor(merged_table$cough)
merged_table$rales <- as.factor(merged_table$rales)


#Imputation

imputed_data <- mice(merged_table, m=5, method = "rf")
summary(imputed_data)
imputed_merged_table <- complete(imputed_data, 1)


#KNN

impute <- kNN(merged_table, variable = c("orthopnea", "oedema", "cough", "rales"), k = 10)


#PCA

install.packages('caret')
library(caret)
merged_table <- na.exclude(merged_table)

smp_size <- floor(0.8 * nrow(merged_table))
set.seed(123)
train_ind <- sample(seq_len(nrow(merged_table)), size = smp_size)
train <- merged_table[train_ind, ]
test <- merged_table[-train_ind, ]

preProc <- preProcess(train[3],method="pca",pcaComp=2)
trainPCA <- predict(preProc, train[3])
model <- train(train$label ~ .,method="glm",data=trainPCA)
testPCA <- predict(preProc,test[3])
predictions_test<-predict(model,testPCA)















