#Load data

#datafolder <- "C:/Users/yyase/Downloads/Core Project Data/"
#datafolder <- "C:/Users/Punkt/Downloads/Core Project Data/"
#datafolder <- "C:/Users/sradu/OneDrive/Documenten/year 3/The core of Biomdical Sciences/Project R/Code Project R/"
datafolder <- "C:/UM/BBS3004 The Core of Biomedical Sciences/Data/"

deaths <- read.csv(paste0(datafolder, "BBS3004_deaths.csv"), header = TRUE)
demo <- read.csv(paste0(datafolder, "BBS3004_demographics.csv"), header = TRUE)
hosp <- read.csv(paste0(datafolder, "BBS3004_hospitalisations.csv"), header = TRUE)
lab <- read.csv(paste0(datafolder, "BBS3004_labvalues.csv"), header = TRUE)
visits <- read.csv(paste0(datafolder, "BBS3004_visits.csv"), header = TRUE)

#Invert the table of labvalues #pivot_wider
.libPaths("C:/Users/Punkt/Downloads/R/RStudio")

library(magrittr)
library(tidyr)

lab_inv <- lab %>%
  pivot_wider(names_from = biomarker, values_from = value)


#Create a common table for visits, labs and demo

merged_table <- data.frame()
merged_table <- merge(x = visits, y = lab_inv)
merged_table <- merge(x = merged_table, y = demo, by = "id")


#hospitalization within 60

hosp <- hosp[hosp$reason == 'Worsening CHF', ]

merged_table$label_hosp <- 0
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
      }
    }
  }
}

      
#Death within 60 days

merged_table$label_death <- 0
for (a in 1:nrow(merged_table)) {
  for (c in 1:nrow(deaths)) {
    if (merged_table$id[a] == deaths$id[c]) {
      if (!is.na(deaths$date_of_death[c]) & difftime((as.Date(deaths$date_of_death[c])),
                     (as.Date(merged_table$date[a])),
                     units = "days") <= 60 &
            difftime((as.Date(deaths$date_of_death[c])),
                     (as.Date(merged_table$date[a])), 
                     units = "days") > 0) {
          merged_table$label_death[a] <- '1'
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
#outliers in NT.BNP --> not normally distributed

outliers <- function(x) {
  
  Q1 <- quantile(x, .25)
  Q3 <- quantile(x, .75)
  IQR <- IQR(x)
  
  outliers <- subset(merged_table, x < (Q1 - 1.5*IQR) | x > (Q3 + 1.5*IQR))
  
  return(outliers)
}

A <- outliers(merged_table$`NT-BNP`)
hist(log(merged_table$`NT-BNP`), 
     main = 'Histogram of NT-BNP distribution',
     xlab = 'log(NT-BNP)',
     col = 'green',
     border = 'blue',
     las = 1)

B <- outliers(merged_table$`CRP sensitive`)
hist(log(merged_table$`CRP sensitive`),
     main = 'Histogram of CRP sensitive distribution',
     xlab = 'log(CRP sensitive)',
     col = 'blue',
     border = 'green',
     las = 1)

C <- outliers(merged_table$`IL-6`)
hist(log(merged_table$`IL-6`),
     main = 'Histogram of Il-6 distribution',
     xlab = 'log(IL-6)',
     col = 'blue',
     border = 'green',
     las = 1)

D <- outliers(merged_table$`GFR`)
hist(log(merged_table$`GFR`),
     main = 'Histogram of GFR distribution',
     xlab = 'log(GFR)',
     col = 'blue',
     border = 'green',
     las = 1)

E <- outliers(merged_table$`Cystatin C`)
hist(log(merged_table$`Cystatin C`),
     main = 'Histogram of Cystatin C distribution',
     xlab = 'log(Cystatin C)',
     col = 'blue',
     border = 'green',
     las = 1)

F <- outliers(merged_table$`age`)
hist(log(merged_table$`age`),
     main = 'Histogram of Age distribution',
     xlab = 'log(age)',
     col = 'blue',
     border = 'green',
     las = 1)

#Update age of the patient
#Date after diagnosis -> days after diagnosis

merged_table$days_after_diagnosis<-difftime(as.Date(merged_table$date), as.Date(merged_table$date_of_diagnosis), units = "days")
merged_table$age <- as.integer(merged_table$age + (merged_table$days_after_diagnosis/365))
merged_table$days_after_diagnosis <- as.integer(merged_table$days_after_diagnosis)
merged_table = merged_table[,!grepl("^date_of_diagnosis",names(merged_table))]


#Imputation of missing values - MICE

#install.packages("mice")
#install.packages("VIM")
#install.packages("Rcpp")
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

#| rename columns because of the spaces
names(merged_table)[7] <- 'B1'
names(merged_table)[8] <- 'B2'
names(merged_table)[9] <- 'B3'
names(merged_table)[10] <- 'B4'
names(merged_table)[11] <- 'B5'

#| remove columns labels and make temporary merge file for imputation
temp_merged <- merged_table[-15]

#Imputation 
imputed_data <- mice(temp_merged, m=5, method = "rf")
summary(imputed_data)
imputed_merged_table <- complete(imputed_data, 1)
#Here I used random forest, maxit and m is on default and I chose model 1. However when we have an actual model we can play around with these and see what results in the best model.

#| Remove id and date of visit
imputed_merged_table <- imputed_merged_table[-c(1:2)]

#PCA & Split

library(caret)
preProc <- preProcess(imputed_merged_table,method="pca",pcaComp=3)
trainPCA <- predict(preProc, imputed_merged_table)

#| No PCA, better AUC spec = 0
#trainPCA <- imputed_merged_table

set.seed(2308)

trainPCA$label <- as.factor(merged_table$label)
intrain <- createDataPartition(y = trainPCA$label, p= 0.8, list = FALSE)
train <- trainPCA[intrain,]
test <- trainPCA[-intrain,]

library(dplyr)
train <- train  %>% 
  mutate(label = factor(label, 
                        labels = make.names(levels(label))))

train_control <- trainControl(method="cv", number=5, classProbs = TRUE, summaryFunction = twoClassSummary)

test <- test  %>% 
  mutate(label = factor(label, 
                        labels = make.names(levels(label))))


#SVM

svm_Linear <- train(label ~., data = train, method = "svmLinear",
                    trControl=train_control,
                    metric = "ROC")
svm_Linear

plot(svm_Linear)

test_pred <- predict(svm_Linear, newdata = test)


#RF

RF <- train(label ~., data = train, method = "rf",
                    trControl=train_control,
                    metric = "ROC")

RF

plot(RF)

test_pred <- predict(RF, newdata = test)


#CART

#install.packages("rpart")
library(rpart)
CART <- train(label ~., data = train, method = "rpart",
            trControl=train_control,
            metric = "ROC")

CART

plot(CART)

test_pred <- predict(CART, newdata = test)


#AUC

#install.packages("pROC")
library(pROC)
roc_obj <- roc(test$label, as.integer(test_pred))
plot(roc_obj)
auc(roc_obj)

#descriptives

colours = c(rep("pink",1), rep("blue", 1))
boxplot(age~gender, data = merged_table, main = 'Boxplot of the average age per gender',
        xlab = 'Gender',
        ylab = 'Age',
        col= colours,
        las = 1)

