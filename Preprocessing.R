#Load data

datafolder <- "C:/Users/yyase/Downloads/Core Project Data/"
#datafolder <- "C:/Users/Punkt/Downloads/Core Project Data/"
#datafolder <- "C:/Users/sradu/OneDrive/Documenten/year 3/The core of Biomdical Sciences/Project R/Code Project R/"
datafolder <- "C:/Users/Punkt/Downloads/Core Project Data/"
#datafolder <- "C:/Users/sradu/OneDrive/Documenten/year 3/The core of Biomdical Sciences/Project R/Code Project R/"
#datafolder <- "C:/UM/BBS3004 The Core of Biomedical Sciences/Data/"

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


#Update age of the patient
#Date after diagnosis -> days after diagnosis

merged_table$days_after_diagnosis<-difftime(as.Date(merged_table$date), as.Date(merged_table$date_of_diagnosis), units = "days")
merged_table$age <- as.integer(merged_table$age + (merged_table$days_after_diagnosis/365))
merged_table$days_after_diagnosis <- as.integer(merged_table$days_after_diagnosis)
merged_table = merged_table[,!grepl("^date_of_diagnosis",names(merged_table))]

#Necessary columns as factor, can also be already done in earlier part of code
merged_table$orthopnea <- as.factor(merged_table$orthopnea)
merged_table$oedema <- as.factor(merged_table$oedema)
merged_table$cough <- as.factor(merged_table$cough)
merged_table$rales <- as.factor(merged_table$rales)

#| rename columns because of the spaces
names(merged_table)[7] <- 'NT_BNP'
names(merged_table)[8] <- 'CRP_sensitive'
names(merged_table)[9] <- 'IL_6'
names(merged_table)[10] <- 'GFR'
names(merged_table)[11] <- 'Cystatin_C'

library(caret)
library(dplyr)

temp_merged <- merged_table
temp_merged$label <- as.factor(temp_merged$label)


#split train and test
intrain <- createDataPartition(y = temp_merged$label, p= 0.8, list = FALSE)
train <- temp_merged[intrain,]
test <- temp_merged[-intrain,]

#only remove id and date of visit, if we do want imputation
train_temp<- train [-c(1:2, 15)]

#imputation
library("mice")
imputed_data <- mice(train_temp, m=5, method = "rf")
summary(imputed_data)
train_temp <- complete(imputed_data, "long")

#add labels back???
label <- rep(train$label, 5)
train <- cbind(train_temp, label)
train$label <- as.factor(train$label)
train<- train[-c(1:2)]

test_temp<- test[-c(1:2, 15)]

#imputation
imputed_data <- mice(test_temp, m=5, method = "rf")
summary(imputed_data)
test_temp <- complete(imputed_data, "long")

#add labels back???
label <- rep(test$label, 5)
test <- cbind(test_temp, label)
test$label <- as.factor(test$label)
test<- test[-c(1:2)]

train <- train  %>% 
  mutate(label = factor(label, 
                        labels = make.names(levels(label))))

train_control <- trainControl(method = "adaptive_cv",
                              number = 5,
                              repeats = 5,
                              adaptive = list (min = 2,
                                               alpha = 0.05,
                                               method ="gls",
                                               complete = FALSE),
                              classProbs = TRUE,
                              summaryFunction = twoClassSummary,
                              verboseIter = TRUE)
                              

test <- test  %>% 
  mutate(label = factor(label, 
                        labels = make.names(levels(label))))


#Classifiers

train_model <- function(x) {
  model <- train(label ~., data = train, method = x,
                 trControl=train_control,
                 metric = "ROC",
                 preProcess = c("center","scale"))
  return(model)
}

predict_model <- function(x) {
  test_pred <- predict(x, newdata = test, type = "prob")
  return(test_pred)
}

library(classifierplots)
calibration <- function(x) {
  calibration_plot(ifelse(test$label == 'X1', 1,0), x$X1) 
}

ROC_AUC <- function (x) {
  classifierplots::roc_plot(ifelse(test$label == 'X1', 1,0), x$X1)
}

#|Logistic Regression
LR<-train_model("glm")
test_LR<-predict_model(LR)
calibration(test_LR)
ROC_AUC(test_LR)

#|SVM
svmLinear<-train_model("svmLinear")
test_svmLinear<-predict_model(svmLinear)
calibration(test_svmLinear)
ROC_AUC(test_svmLinear)

#|RF
RF<-train_model("rf")
test_RF<-predict_model(RF)
calibration(test_RF)
ROC_AUC(test_RF)
plot(varImp(RF, scale = TRUE))
plot(RF)

#|CART
library(rpart)
CART<-train_model("rpart")
test_CART<-predict_model(CART)
calibration(test_CART)
ROC_AUC(test_CART)

#|kNN
kNN<-train_model("knn")
test_kNN<-predict_model(kNN)
calibration(test_kNN)
ROC_AUC(test_kNN)

#|ANN
ANN<-train_model("nnet")
test_ANN<-predict_model(ANN)
calibration(test_ANN)
ROC_AUC(test_ANN)

#|NB
NB<-train_model("nb")
test_NB<-predict_model(NB)
calibration(test_NB)
ROC_AUC(test_NB)









