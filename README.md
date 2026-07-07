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


or_table <- tidy(model, conf.int = TRUE, exponentiate = TRUE) 

or_table


