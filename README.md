#1.Graphically drawing the lines;
x<-seq(0,400,1)
yc<-
yd<-
yl<-
plot(x, yc, type="n",xlab="x",ylab="y",
     xlim=c(0,400), ylim=c(0,400))
lines(x,yc,col="red")
lines(x,yd,col="blue")
lines(x,yl,col="green")

The corner points are: (0, 120), (160, 0), (120, 40).
Calculating the objective value for all corner points;

df<-data.frame(x=c(0,160,120),y=c(120,0,40))
df$z<-100*df$x+120*df$y
opt<-which.max(df$z)
xy<-df[opt,]
xy

#2 simplex method

library(lpSolve)
obj<-c()
const.mat<-matrix(c(),nrow = 3,byrow = TRUE)
const.dir<-c("<=","<=","<=")
const.rhs<-c()
sol<-lp("min",obj,const.mat,const.dir,const.rhs)
cat("Z=",sol$objval,"\n")
cat("(x1, x2,x3):",sol$sol[1:3],"\n")


# dual simplex method (To solve the given LPP by dual simplex method we have to maximize the
objective function and have to change the constraint into less than or equal.
)
library(lpSolve)
obj<-c()
const.mat<-matrix(c(),nrow = 2,byrow = TRUE)
const.dir<-c("<=","<=")
const.rhs<-c(-2,-4)
sol<-lp("max",obj,const.mat,const.dir,const.rhs)
cat("Z=",sol$objval,"\n")
cat("(x1, x2):",sol$sol[1:2],"\n")

#Big-M method(To solve the given LPP by Big-M method, we write the given LPP introducing slack, surplus and artificial variables)

library(lpSolve)
M<-1e10
obj<-c(3,2,0,0,0,M,M)
const.mat<-matrix(c(2,1,0,-1,0,1,0,-3,2,1,0,0,0,0,1,1,0,0,-1,0,1),nrow = 3,byrow = TRUE)
const.dir<-c("=","=","=")
const.rhs<-c(10,6,6)
sol<-lp("min",obj,const.mat,const.dir,const.rhs)
cat("Z=",sol$objval,"\n")
cat("Opt Values (x1, x2):",sol$sol[1:2],"\n")

# two-phase method()
#p1
library(lpSolve)
obj1<-c(0,0,0,0,1,1)
const1<-matrix(c(2,1,1,1,3,1,1,0,0,0,-1,0,0,1,0,0,0,1),nrow=3)
dir1<-c("=","=","=")
rhs1<-c(16,36,10)
phase1<-lp(
  direction="min",
  objective.in=obj1,
  const.mat=const1,
  const.dir=dir1,
  const.rhs=rhs1)
cat("the optimal value =",phase1$objval,"\n")

cat("The values of the artificial variables are:","\n",
    "a1=",phase1$solution[5]," and ","a2=",phase1$solution[6],"\n")
#p2
if (phase1$objval == 0) {
  obj.phase2 <- c(2, 3, 0, 0, 0, 0)
  phase2 <- lp("min", obj.phase2, const.mat, const.dir, const.rhs)
  cat("\nPhase 2 (original objective):\n")
  cat("Optimal value of Z:", phase2$objval, "\n")
  cat("Values of x1, x2:", phase2$solution[1:2], "\n")
} else {
  cat("Problem is infeasible (Phase 1 > 0)\n")
}


#transportation

library(lpSolve)
m<-matrix(c(),nrow=3,byrow=T)
a<-c("<","<","<")
b<-c()
c<-c(">",">",">",">")
d<-c()
sol<-lp.transport(m,"min",a,b,c,d)
sol
sol$solution

#NCWR and VAM
# Define the cost matrix
costs <- matrix(c(4, 8, 8,
                  16, 24, 16,
                  8, 16, 24),
                nrow = 3, byrow = TRUE)
supply <- c(76, 82, 77)
demand <- c(72, 102, 41)
# Function for North-West Corner Rule
north_west_corner <- function(supply, demand) {
  m <- length(supply)
  n <- length(demand)
  allocation <- matrix(0, m, n)
  
  i <- 1
  j <- 1
  while (i <= m & j <= n) {
    alloc <- min(supply[i], demand[j])
    allocation[i, j] <- alloc
    supply[i] <- supply[i] - alloc
    demand[j] <- demand[j] - alloc
    
if (supply[i] == 0) {
      i <- i + 1
    } else {
      j <- j + 1
    }
  }
  return(allocation)
}
# Run NWCR
nwcr_result <- north_west_corner(supply, demand)
cat("North-West Corner Rule Allocation:\n")
print(nwcr_result)

# Function for Vogel’s Approximation Method (VAM)
vogel_approximation <- function(costs, supply, demand) {
  m <- nrow(costs)
  n <- ncol(costs)
  allocation <- matrix(0, m, n)
  supply_left <- supply
  demand_left <- demand
  cost_matrix <- costs
  
  while (any(supply_left > 0) && any(demand_left > 0)) {
    penalty <- rep(0, m + n)
    
    # Calculate row penalties
for (i in 1:m) {
      if (supply_left[i] > 0) {
        row <- cost_matrix[i, demand_left > 0]
        if (length(row) >= 2) {
          penalty[i] <- sort(row)[2] - sort(row)[1]
        } else if (length(row) == 1) {
          penalty[i] <- row
        }
      }
    }
    
    # Calculate column penalties
for (j in 1:n) {
      if (demand_left[j] > 0) {
        col <- cost_matrix[supply_left > 0, j]
        if (length(col) >= 2) {
          penalty[m + j] <- sort(col)[2] - sort(col)[1]
        } else if (length(col) == 1) {
          penalty[m + j] <- col
        }
      }
    }
    
    # Find max penalty
idx <- which.max(penalty)
    
if (idx <= m) {
      # Row penalty selected
i <- idx
      available_cols <- which(demand_left > 0)
      j <- available_cols[which.min(cost_matrix[i, available_cols])]
    } else {
      # Column penalty selected
      j <- idx - m
      available_rows <- which(supply_left > 0)
      i <- available_rows[which.min(cost_matrix[available_rows, j])]
    }
    
alloc <- min(supply_left[i], demand_left[j])
    allocation[i, j] <- alloc
    supply_left[i] <- supply_left[i] - alloc
    demand_left[j] <- demand_left[j] - alloc
  }
  return(allocation)
}
# Run VAM
vam_result <- vogel_approximation(costs, supply, demand)
cat("\nVogel’s Approximation Method Allocation:\n")
print(vam_result)

#####bio#####

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

#############
# Example-07:Kaplen Meier Estimate
# ----------------------------------------------------------------
library(readr)
SP <- read_csv("C:/Users/User/Desktop/Schizophrenia.csv")

str(SP)

SP <- within(SP,{
  cen <- factor(Censor, labels = c("Censor","Death"))
  marital <- factor(Marital, labels = c("Single","Married","Alone again"))
  sex <- factor(Gender, labels = c("Male","Female"))
})

str(SP)

# Number of censored observations and events for all patients: 
table(SP$cen)
prop.table(table(SP$cen))


# Number of censored observations and events for all patients(by marital status): 
table(SP$marital,SP$cen)
prop.table(table(SP$marital,SP$cen))

# Number of censored observations and events for all patients(by gender): 
table(SP$sex,SP$cen)
prop.table(table(SP$sex,SP$cen))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kaplan-Meier survival curve for all patients
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(survival)
names(SP)
sp_survival <- survfit(Surv(Time, Censor)~1, data = SP, type = "kaplan-meier", conf.type = "log-log")

library(survminer)
ggsurvplot(sp_survival,
           data = SP,
           risk.table = TRUE,
           pval = FALSE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The median survival time together with a 95% CI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sp_survival


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The survival time at 6 months, 18 months, and 36 months together with a 95% CI
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary(sp_survival, time = c(6*30,18*30,36*30))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The median survival time (with 95% CI) by marital status
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(SP)
sp_survival_mar <- survfit(Surv(Time,Censor)~marital, data = SP, type = "kaplan-meier", conf.type = "log-log")

sp_survival_mar 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kaplan-Meier estimate of the survival curve by marital status
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggsurvplot(sp_survival_mar,
           data = SP,
           risk.table = TRUE,
           pval = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The median survival time (with 95% CI) by gender
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(SP)
sp_survival_sex <- survfit(Surv(Time,Censor)~sex, data = SP, type = "kaplan-meier", conf.type = "log-log")

sp_survival_sex 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Kaplan-Meier estimate of the survival curve by gender
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggsurvplot(sp_survival_sex,
           data = SP,
           risk.table = TRUE,
           pval = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Comparison with Nelson-Aalen estimate at timepoints 90, 180, 270, 365 days
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


sp_km <-survfit(Surv(Time,Censor)~1, data = SP, type = "kaplan-meier", conf.type = "log-log")      # Kaplan-Meier
sp_na <-survfit(Surv(Time,Censor)~1, data = SP, type = "fh", conf.type = "log-log")                # Nelson Aalen / Flemming Harington

sum_km <- summary(sp_km, times = c(90, 180, 270, 365))
sum_na <- summary(sp_na, times = c(90, 180, 270, 365))

data.frame(
  Time = sum_km$time,
  At_risk = sum_km$n.risk,
  Events = sum_km$n.event,
  Censor = sum_km$n.censor,
  Skm =sum_km$surv,
  km_lower = sum_km$lower,
  km_upper = sum_km$upper,
  Sna = sum_na$surv,
  na_lower = sum_na$lower,
  na_upper = sum_na$upper
)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Confidence interval types: (simple, log, clog-log)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sp_surv_plain <- survfit(Surv(Time,Censor)~1, data = SP, type = "kaplan-meier", conf.type = "plain")   # simple
summary(sp_surv_plain)

sp_surv_log <- survfit(Surv(Time,Censor)~1, data = SP, type = "kaplan-meier", conf.type = "log")   # log
summary(sp_surv_log)

sp_surv <- survfit(Surv(Time,Censor)~1, data = SP, type = "kaplan-meier", conf.type = "log-log")   # clog-log
summary(sp_surv)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Taking variance from std error
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sp_surv <- survfit(Surv(Time,Censor)~1, data = SP, type = "kaplan-meier", conf.type = "log-log")
sum <- summary(sp_surv,  times = c(90, 180, 270, 365))

data.frame(
  time = sum$time,
  nj = sum$n.risk,
  dj = sum$n.event,
  cj = sum$n.censor,
  S = sum$surv,
  Var = sum$std.err^2
)


#example-03
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Type I censoring 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

lung <- survival::lung  # Importing dataset

lung$status2 <- ifelse(lung$status==2,1,0)
n <- nrow(lung)

names(lung)

L <- 365

ti <- pmin(lung$time,L)
si <- as.numeric(lung$status2 == 1 & lung$time <= L)

r <- sum(si)
T <- sum(ti)

cat("\nUnder type I censoring\n Censoring time:",L,"\n",
    "Number of events:",r,"|Number of censor:",n-r)

lam_hat <- r/T
mu_hat <-  1/lam_hat

lam_hat
mu_hat

sd_lam <- lam_hat/sqrt(r)
round(lam_hat+c(-1,1)*1.96*sd_lam,4)

sd_mu <- mu_hat/sqrt(r)
round(mu_hat+c(-1,1)*1.96*sd_mu,4)

qchisq(0.025,2*r)/(2*T)
qchisq(0.975,2*r)/(2*T)



# ~~~~~~~~~~~~~~
# Example-4:estimation of standard error
# ~~~~~~~~~~~~~~
time1 <- c(1,1,3,5,5,14,17,17,23,23)
time2 <- c(2,2,6,9,9,10,11,12,13,13,13,14,17,18,19,21,21,23,24,24) 

event <- c(rep(1,10),rep(0,20))
time <- c(time1, time2)

df <- data.frame(
  Event = event,
  Time = time
)

df

library(survival)
mod <- survfit(Surv(Time, Event)~1, data = df, type = "kaplan-meier", conf.type = "log-log")

mod      #median survival year/median remission 95% CI
summary(mod)

#probability that participants survive 16 and 21 years
summary(mod, times = c(16,21))


#survival curve
library(survminer)
ggsurvplot(mod,
           data = df,
           risk.table = TRUE,
           pval = FALSE,
           xlab = "days to follow-up",
           conf.int = FALSE)

s16 <- summary(mod, time = 14)


#if null hypothesis is given by s(t16)=0.6
s16e <- s16$surv
s16e
s16sd <- s16$std.err
s16sd

zcal <- (s16e-0.6)/s16sd  #here 0.6 from ques
zcal
ztab <- qnorm(0.975)
ztab
srt(time1)

#test for equality of the survival function in two group
time1$group<-"time1"
time2$group<-"time2"
time<-rbind(time1,time2)
library(survival)
survdiff(Surv(time,event)~group,data=time)



# ~~~~~~~~~~~~~~~~
# Example 5 -----
# ~~~~~~~~~~~~~~~~

time1 <- c(1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,6,6,9,9,9,10,13,16)
time2 <- c(2.1,3.1,4,7,7,8,8,9,9,9,11,24,24)

time <- c(time1, time2)
event <- c(rep(1,25),rep(0,13))

df <- data.frame(time, event)

df
str(df)

mod <- survfit(Surv(time,event)~1,df,type="kaplan-meier",conf.type="log-log")
mod
summary(mod)

summary(mod)$n.censor


# ~~~~~~~~~~~~~~~~
# Example 6 -----
# ~~~~~~~~~~~~~~~~

time <- c(6,6,6,7,10,13,16,22,23,6,9,10,11,17,19,20,25,32,32,34,35)
event <- c(rep(1,9),rep(0,12))

df1 <- data.frame(time,event)

time2 <- c(1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)
event2 <- rep(1,21)

df2 <- data.frame(time2,event2)

library(survival)
library(survminer)
s1 <- survfit(Surv(time,event)~1, 
              df1, type = "kaplan-meier", conf.type = "log-log") 

s2 <- survfit(Surv(time2,event2)~1, 
              df2, type = "kaplan-meier", conf.type = "log-log") 

#Survival Estimates for group 1 and group 2 with variance:
r1 <- summary(s1)

r2 <- summary(s2)


data.frame(r1$time, r1$n.risk, r1$n.event, r1$n.censor, 
           round(r1$surv,4), variance = round(r1$std.err^2,4))

data.frame(r2$time, r2$n.risk, r2$n.event, r2$n.censor, 
           round(r2$surv,4), variance = round(r2$std.err^2,4))

#Merging to groups to plot two curve in the same graph:

group <- c(rep(1,21),rep(2,21))
df <- data.frame(Time = c(time,time2),Event =c(event,event2),group)

s <- survfit(Surv(Time,Event)~group, 
             df, type = "kaplan-meier", conf.type = "log-log")


ggsurvplot(s,
           data = df,
           risk.table = TRUE,
           pval = TRUE,
           conf.int = FALSE,
           surv.median.line = "hv",
           xlab = "Remision times(in weeks)")

# 95% CI of S(t) for group 1 at time = 16:

summary(s1,times=16)

# 95% confidence interval for the median remission time for group-1:

s1   # Answer = (13, Not estimable)


# Test the null hypothesis:
# H0: s(t=16) = 0.9

# Test statistic, z = s_hat(t=16)-s(t=16) / se(s_hat(t=16))
# Not that under H0, s(t=16) = 0.9

s16 <- summary(s1,times = 16)
shat_16 <- s16$surv
se <- s16$std.err

z_cal <- (shat_16-0.9)/se

z_tab <- qnorm(0.975)   #Alpha = 5%

z_cal
z_tab

abs(z_cal)>abs(z_tab) # TRUE, reject the null hypothesis.

# Equality of the survival functions in the two groups

lrt <- survdiff(Surv(Time,Event)~group,df)

lrt  # p value <0.001, so there exist a difference in survival between two groups.






# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Logistic Distribution ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# 1) Importing data set:
library(readr)
HT <- read_csv("D:/Study Files/8th Sem/Bio lab/Hypertension_Lab.csv")
sum(is.na(HT$hypert))

# 2) Re-coding variables:
str(HT)

HT <- within(HT,{
  bmi_cat <- factor(ifelse(bmi<25, "Normal", 
                           ifelse(bmi<30, "Overweight", "Obese")),
                    levels = c("Normal","Overweight", "Obese"))
  sex <- factor(Gender, labels = c("Male","Female"))
  ht <- factor(hypert, labels = c("No", "Yes"))
})

str(HT)


# Table and proportion table of gender:
t <- table(HT$sex, HT$ht)
t

round(prop.table(t),4)

# Risk of hypertension for male:
t

risk <- t[1,2]/(t[1,1]+t[1,2])
risk 

# 95% CI for the risk:

risk+c(-1,1)*qnorm(0.975)*sqrt(risk*(1-risk)/(t[1,1]+t[1,2]))


# Association between the hypertension and bmi category:
t2 <- table(HT$bmi_cat, HT$ht)

chisq.test(t2)   # Since p-value = 1.011e-14 <0.001, there is a significant association.

# Highest probability of hypertension

round(prop.table(t2,2),4)*100 # 2 means column wise percentage, so the category which has higher percentage for ht has the highest probability. 


# Fit a logistic regression model to predict hypertension, using bmi, sleep, gender:

mod <- glm(ht ~ bmi + sex + sleep, data = HT, family = binomial(link = "logit"))
summary(mod)

library(broom)
tidy(mod, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)

# Deviance of the model:
cat("\nNull Devience:",mod$null.deviance,
    "\nResidual Devience:", mod$deviance)

# Calculate the sensitivity, recall, specificity, PPV, precision, NPV, F1 score and accuracy level of the test and interpret:
pp <- predict(mod, type = "response")
pv <- factor(ifelse(pp>0.5, 1, 0), labels = c("No", "Yes"))

ct <- table(Actual = HT$ht,
            Predicted = pv)
ct

TN <- ct[1,1]
FP <- ct[1,2]
FN <- ct[2,1]
TP <- ct[2,2]

calc_mat <- function(TP,TN,FP,FN){
  Total <- TP+TN+FP+FN
  Sensitivity <- TP/(TP+FN)
  Specificity <- TN/(FP+TN)
  PPV <- TP/(TP+FP)
  NPV <- TN/(TN+FN)
  Accuracy <- (TP+TN)/Total
  
  Recall <- Sensitivity
  Precision <- PPV
  F1 <- 2*Recall*Precision/(Recall+Precision)
  
  t(data.frame(Total,
               Sensitivity,
               Specificity,
               PPV,
               NPV,
               Accuracy,
               Recall,
               Precision,
               F1))
}

calc_mat(TP,TN,FP,FN)

# 24.79% of the actual cases of positive hypertension is predicted by the model.
# 93.20% of the actual cases of negative hypertension is predicted by the model.
# If the model predicts an individual having hypertension then there is a 64.52% chance of actually having hypertension.
# If the model predicts an individual not having hypertension then there is a 71.29% chance of actually not having hypertension.
# The model correctly predicts 70.43% of the actual cases.
# Since F1 score = 0.35 is means the model performs poorly.


# Multicollinearity by variance inflation factor
library(car)
vif(mod)

# VIF indicator: vif = 1 ; No
#               vif  >=5  ; Severe 


# Hosmer-Lemeshow goodness of fit:

# H0: The model is a good fit
# H1: The model is not a good fit.

library(ResourceSelection)
hoslem.test(x = as.numeric(HT$ht)-1,
            y = fitted(mod),
            g=10)

# Since p value>0.05, hence it's a good fit.


# McFadden's pseudo-R2:
mcfad <- 1-mod$deviance/mod$null.deviance
mcfad

cat("\nThe model explains",round(mcfad,4)*100,"percent of the variation in the response variable.
    \nIt's a poor fit.")


# ROC Curve/AUC:
library(pROC)

roc_obj <- roc(response = HT$ht,
               predict = fitted(mod),
               levels = c("No","Yes"))
auc(roc_obj)

plot(roc_obj, print.auc=TRUE)

cat("\nIf we randomly pick one positive case and one negative case then the model will assign higher probablity to 
    the positive case about",round(auc(roc_obj),4)*100,"percent of the time.")



# Average Marginal Effects:
library(margins)
ame <- margins(mod)
summary(ame)


or_table <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 

or_table


