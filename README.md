
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse,        
  marginaleffects,  
  car,            
  ResourceSelection,
  pROC,             
  broom,          
  ggplot2
)

setwd("F:/Statistics semester wise book & sheet/8th semester/STA 4205/BioStat/")
HT <- read.csv("Hypertension_Lab.csv",header = T, sep = ",")
names(HT)


HT$bmi_cat <- ifelse(HT$bmi < 25, "Normal",
                        ifelse(HT$bmi < 30, "Overweight","Obese"))
      # Convert to factor
HT$bmi_cat <- factor(HT$bmi_cat, 
                        levels = c("Normal", "Overweight", "Obese"))
HT <- within(HT, {
  sex <- factor(Gender, labels = c("Male", "Female"))
  ht <- factor(hypert, labels = c("No","Yes"))
  })
str(HT)

prop.table(table(HT$sex))


 
HT |>
  group_by(ht) |>
  summarise(
    n          = n(),
    mean_bmi   = mean(bmi),
    sd_bmi     = sd(bmi),
    mean_sleep = mean(sleep),
    sd_sleep   = sd(sleep),
    ) |>
  print()

table(HT$ht,HT$sex)

model <- glm(
  ht ~ bmi + sleep + sex,
  data   = HT,
  family = binomial(link = "logit")
)

summary(model)


tidy(model, conf.int = TRUE, conf.level = 0.95) |> print()




predicted_prob <- predict(model, type="response")

predicted_class<-ifelse(predicted_prob > 0.5, 1, 0)

classification_table<- table(Actual = HT$ht,
                            Predicted = predicted_class)
print(classification_table)

TN <- classification_table[1, 1]  #  True Negatives (452)
FP <- classification_table[1, 2]  #  False Positives (33)
FN <- classification_table[2, 1]  #  False Negatives (182)
TP <- classification_table[2, 2]   #  True Positives (60)

calculate_metrics<- function (TP, TN, FP, FN) {
  metrics <- list(
  TP <- TP,
  TN <- TN,
  FP <- FP,
  FN <- FN,
  Total <- TP + TN + FP + FN,
  Accuracy <- (TP + TN)/(TP + TN + FP + FN),
  Sensitivity <- TP/(TP + FN),
  Recall <- TP/(TP+ FN),
  Specificity <- TN/(TN +FP),
  PPV <- TP/(TP + FP),
  Precision <- TP/(TP + FP),
  NPV <- TN/(TN + FN),
  F1_Score <- 2 * (TP/(TP+FP)) * (TP/(TP+FN)) / ((TP/(TP+FP)) + (TP/(TP+FN)))
  )
}

results <- calculate_metrics (TP, TN, FP, FN)
print(results)


cat("\nVIF:\n")
print(vif(model))


hl <- hoslem.test(
  x = as.numeric(HT$ht) - 1,
  y = fitted(model),
  g = 10
)

cat("\nHosmer-Lemeshow test:\n")
print(hl)

cat("\nNull deviance   :", model$null.deviance,  "df =", model$df.null)
cat("\nResidual deviance:", model$deviance,       "df =", model$df.residual)
cat("\nAIC             :", AIC(model), "\n")

mcfadden_r2 <- 1 - (model$deviance / model$null.deviance)
cat("\nMcFadden R²:", round(mcfadden_r2, 4), "\n")

library(pROC)
roc_obj <- roc(
  response  = HT$ht,
  predictor = fitted(model),
  levels    = c("No", "Yes")
)
cat("AUC:", round(auc(roc_obj), 4), "\n")

plot(roc_obj,
     col  = "#4F46E5",
     lwd  = 2,
     main = paste0("ROC Curve  (AUC = ", round(auc(roc_obj), 3), ")"),
     print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "gray60")



ame <- avg_slopes(model)
print(ame)
names(ame)

cat("\n--- AME Interpretation ---\n")
cat("Each 1-unit increase in BMI changes P(hypertension) on average by",
    round(ame$estimate[ame$term == "bmi"], 4), "\n")
cat("Being Female (vs male) changes P(hypertension) on average by",
    round(ame$estimate[ame$term == "genderFemale"], 4), "\n")
cat("Each 1-hour increase in sleep changes P(hypertension) on average by",
    round(ame$estimate[ame$term == "sleep"], 4), "\n")

mem <- slopes(model, newdata = datagrid())
print(mem)


or_table <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 

or_table


