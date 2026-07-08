#1a
x<-seq(0,400,1)
yc<-(34000-100*x)/100
yd<-(36000-200*x)/300
yl<-
plot(x, yc, type="n",xlab="x",ylab="y",
     xlim=c(0,400), ylim=c(0,400))
lines(x,yc,col="red")
lines(x,yd,col="blue")
lines(x,yl,col="green")
df<-data.frame(x=c(0,160,120,0),y=c(0,0,40,120))
df$z<-100*df$x+120*df$y
opt<-which.max(df$z)
xy<-df[opt,]
xy
#1b
library(lpSolve)
obj<-c(1,1,-4)
const.mat<-matrix(c(1,1,-1,1,1,1,2,-1,1),nrow=3)
const.dir<-c("<=","<=","<=")
const.rhs<-c(8,2,4)
sol<-lp("min",obj,const.mat,const.dir,const.rhs)
cat("Z=",sol$objval,"\n")
cat("(x1, x2,x3):",sol$sol[1:3],"\n")
#2a
library(lpSolve)
M<-1e10
obj<-c(3,2,0,0,0,M,M)
const.mat<-matrix(c(2,1,0,-1,0,1,0,-3,2,1,0,0,0,0,1,1,0,0,-1,0,1),nrow = 3,byrow = TRUE)
const.dir<-c("=","=","=")
const.rhs<-c(10,6,6)
sol<-lp("min",obj,const.mat,const.dir,const.rhs)
cat("Z=",sol$objval,"\n")
cat("Opt Values (x1, x2):",sol$sol[1:2],"\n")


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
#2b

#NCWR and VAM
cost <- matrix(c(4,8,8,0,16,24,16,0,8,16,24,0),3,4,by=T)
sup <- c(76,82,77)
sum(sup)
dem <- c(72,102,41,20)
sum(dem)
# NWCR
NWCR <- function(s,d){
  a <- matrix(0,length(s),length(d)); i <- j <- 1
  while(i<=length(s) & j<=length(d)){
    x <- min(s[i],d[j]); a[i,j] <- x
    s[i] <- s[i]-x; d[j] <- d[j]-x
    if(s[i]==0) i <- i+1 else j <- j+1
  }
  a
}
NWCR(sup,dem)
VAM <- function(c,s,d){
  m <- nrow(c); n <- ncol(c); a <- matrix(0,m,n)
  while(any(s>0)&any(d>0)){
    p <- c(
      sapply(1:m,\(i) if(s[i]){x<-sort(c[i,d>0]); if(length(x)>1)x[2]-x[1] else x}else 0),
      sapply(1:n,\(j) if(d[j]){x<-sort(c[s>0,j]); if(length(x)>1)x[2]-x[1] else x}else 0)
    )
    k <- which.max(p)
    if(k<=m){i<-k;j<-which(d>0)[which.min(c[i,d>0])]}
    else{j<-k-m;i<-which(s>0)[which.min(c[s>0,j])]}
    x <- min(s[i],d[j]); a[i,j] <- x
    s[i] <- s[i]-x; d[j] <- d[j]-x
  }
  a
}
VAM(cost,sup,dem)

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

#######
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
setwd("")
SP <- read.csv("Schizophrenia.csv",header = T, sep = ",")

str(SP)
view(SP)

SP <- within(SP,{
  cen <- factor(Censor, labels = c("Censor","Death"))
  marital <- factor(Marital, labels = c("Single","Married","Alone again"))
  sex <- factor(Gender, labels = c("Male","Female"))
})

str(SP)
# a Number of censored observations and events for all patients
table(SP$cen)
prop.table(table(SP$cen))

#(by marital status): 
table(SP$marital,SP$cen)
prop.table(table(SP$marital,SP$cen))

# (by gender): 
table(SP$sex,SP$cen)
prop.table(table(SP$sex,SP$cen))
#b km plot
library(survival)
names(SP)
sp_survival <- survfit(Surv(Time, Censor)~1,
                       data = SP, type = "kaplan-meier",
                       conf.type = "log-log")
plot(sp_survival, xlab = "Days of follow-up", ylab="Survival Probability", main= "Overall
survival curve")

#or
library(survminer)
ggsurvplot(sp_survival,
           data = SP,
           risk.table = TRUE,
           pval = FALSE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")

# c median s(t)
print(sp_survival)
# 95% CI
summary(sp_survival, time = c(6*30,18*30,36*30))

#d 95% CI by ms
names(SP)
sp_survival_mar <- survfit(Surv(Time,Censor)~marital, data = SP, type = "kaplan-meier", conf.type = "log-log")

sp_survival_mar 

#plot km by ms
ggsurvplot(sp_survival_mar,
           data = SP,
           risk.table = TRUE,
           pval = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")
cat("\nSince p-value < 0.0001, we conclude that highly significant differences in survival
times among the three marital status groups. Patients who are 'Alone Again' have the poorest survival (median 539 days), while married patients have the best
survival (median 1310 days).") 





# e The median survival time (with 95% CI) by gender
names(SP)
sp_survival_sex <- survfit(Surv(Time,Censor)~sex, data = SP, type = "kaplan-meier", conf.type = "log-log")

sp_survival_sex 
cat("Females have a substantially lower median survival time compared to males,
indicating that females in this cohort have poorer survival outcomes.")

# Kaplan-Meier estimate of the survival curve by gender
ggsurvplot(sp_survival_sex,
           data = SP,
           risk.table = TRUE,
           pval = TRUE,
           surv.median.line = "hv",
           conf.int = FALSE,
           xlab = "Days to follow-up")




cat("The Kaplan-Meier curves show that females have consistently lower survival probabilities than males throughout the follow-up period, with a log-rank p-value <
0.0001 confirming a statistically significant difference")





#f Log-rank test for group comparison
logrank_test <- survdiff(
  Surv(Time, Censor) ~ marital,
  data = SP
)

print(logrank_test)

#4

#HT
setwd("G:/STATISTICS/4rth Year/2nd Semester/lab")
HT <- read.csv("Hypertension_Lab.csv",header = T, sep = ",")

#a
HT$age<-HT$ID/HT$bmi
view(HT)
#c factorize
HT <- within(HT,{
  bmi_cat <- factor(ifelse(bmi<25, "Normal", 
                           ifelse(bmi<30, "Overweight", "Obese")),
                    levels = c("Normal","Overweight", "Obese"))
  sex <- factor(Gender, labels = c("Male","Female"))
  ht <- factor(hypert, labels = c("No", "Yes"))
})

str(HT)


#d Table and proportion table of gender:
t <- table(HT$sex, HT$ht)
t

round(prop.table(t),4)


#f Association between the hypertension and bmi category:
t2 <- table(HT$bmi_cat, HT$ht)

chisq.test(t2)                                                                    # Since p-value = 1.011e-14 <0.001, there is a significant association.

# Highest probability of hypertension

round(prop.table(t2,2),4)
cat("Obese bmi category has the highest probability of hypertension.")

#g Fit a logistic regression model to predict hypertension, using bmi, sleep, gender:

mod <- glm(ht ~ bmi + sex + sleep, data = HT, family = binomial(link = "logit"))
summary(mod)


#h,i
library(broom)
tidy(mod, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95)

Holding sleep duration and gender constant, each one-unit increase in BMI
increases the odds of hypertension by approximately 11%.
After adjusting for BMI and gender, each additional hour of sleep decreases the
odds of hypertension by about 1.4%.


#k Calculate the sensitivity, recall, specificity, PPV, precision, NPV, F1 score and accuracy level of the test and interpret:
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












cat("\nSensitivity (24.79%) / Recall: Only 24.79% of hypertensive individuals were correctly identified.
Specificity (93.20%): 93.20% of non-hypertensive individuals were correctly identified.
PPV / Precision (64.52%): Among those predicted to have hypertension, 64.52% actually had hypertension.
NPV (71.29%): Among those predicted to be non-hypertensive, 71.29% were truly non-hypertensive.
F1 Score (35.82%): Indicates poor overall performance in detecting hypertension.
Accuracy (70.43%): The model correctly classified 70.43% of individuals.")

#p ROC/AUC curve:
library(pROC)

roc_obj <- roc(response = HT$ht,
               predict = fitted(mod),
               levels = c("No","Yes"))
auc(roc_obj)








cat("\nApproximately 68.5% probability that the model assigns a higher predicted
probability of hypertension to a randomly selected hypertensive individual than to
a randomly selected non-hypertensive individual.")

plot(roc_obj, print.auc=TRUE)






cat("\nIf we randomly pick one positive case and one negative case then the model will assign higher probablity to 
    the positive case about",round(auc(roc_obj),4)*100,"percent of the time.")



#q Average Marginal Effects:
library(margins)
ame <- margins(mod)
summary(ame)








cat("\nHolding all other variables constant, each one-unit increase in BMI increases the
probability of hypertension by approximately 2.08 percentage points on average.
After controlling for BMI and sleep, females have an average 3.91 percentage-point
lower probability of hypertension than males.
Holding other variables constant, an additional hour of sleep decreases the
probability of hypertension by approximately 0.29 percentage points on average.")

