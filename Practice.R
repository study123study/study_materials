# Packages you need
if (!requireNamespace("pacman", quietly = TRUE))
  install.packages("pacman")

pacman::p_load(
  tidyverse,         
  marginaleffects,    
  car,               
  ResourceSelection,  
  pROC,              
  broom,              
  ggplot2
)
setwd("G:/STATISTICS/4rth Year/2nd Semester/lab")
HT <- read.csv("Hypertension_Lab.csv",header = T, sep = ",")
names(HT)
view(HT)
sum(is.na(HT))
HT$bmi_cat <- ifelse(HT$bmi < 25, "Normal",
                 ifelse(HT$bmi < 30, "Overweight","Obese"))
HT <- within(HT, {
  bmi_cat <- factor(bmi_cat, levels = c("Normal", "Overweight", "Obese"))
  sex <- factor(Gender, labels = c("Male", "Female"))
  ht <- factor(hypert, labels = c("No","Yes"))
})
prop.table(table(HT$sex))
table(HT$ht, HT$sex)
t<-table(HT$bmi_cat, HT$ht)
t
chisq.test(t)
prop.table(t,2)
model <- glm(ht~bmi+sleep+sex, data=HT, family=binomial(link= "logit"))
library(broom)
tidy(model,exponentiate=T)
cat("\nNull deviance :", model$null.deviance)
cat("\nResidual deviance:", model$deviance)
cat("\nAIC :", AIC(model))
predicted_prob <- predict(model, type="response")
predicted_class<-ifelse(predicted_prob > 0.5, 1, 0)
classification_table<- table(Actual = HT$ht, Predicted = predicted_class)
classification_table
TN <- classification_table[1, 1]
FP <- classification_table[1, 2]
FN <- classification_table[2, 1]
TP <- classification_table[2, 2]
sen<-TP/(TP+FN)
rec<-TP/(TP+FN)
spc<-TN/(TN+FP)
ppv<-TP/(TP+FP)
prec<- TP/(TP + FP)
npv<-TN/(TN+FN)
f1<- 2*(TP/(TP+FP))*(TP/(TP+FN))/((TP/(TP+FP))+(TP/(TP+FN)))
acc<-(TP+TN)/(TP+TN+FP+FN)
cat("\n Sensitivity =",sen,"\n", "Recall =",rec, "\n", "Specificity =",spc, "\n", "PPV
=",ppv,"\n", "Precision =",prec, "\n", "NPV =",npv,"\n", "F1 score =",f1, "\n",
    "Accuracy =",acc,"\n")
library(car)
vif(model)
library(ResourceSelection)
hoslem.test(x = as.numeric(HT$hypert), y = fitted(model), g = 10)
1-(model$deviance/model$null.deviance)
library(pROC)
roc_obj <- roc(response = HT$ht, predictor = fitted(model), levels = c("No",
                                                                       "Yes"))
plot(roc_obj)
auc(roc_obj)
library(marginaleffects)
ame <- avg_slopes(model)
ame
