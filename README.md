#ANALYZE
treat<-factor(rep(paste0("Diet_",1:6),c(6,7,8,6,8,5)))
d1<-c(179,186,200,182,187,197)
d2<-c(169,188,179,189,186,168,182)
d3<-c(182,191,243,214,213,202,180,214)
d4<-c(202,220,226,245,241,225)
d5<-c(210,214,220,192,188,189,190,224)
d6<-c(155,176,159,161,165)
diet<-c(d1,d2,d3,d4,d5,d6)
data<-data.frame(Treatment=treat, Weight=diet)
model<-aov(Weight~Treatment, data=data)
summary(model)
#3&5
con<-matrix(c(0,0,1,0,-1,0,0,0,1,0,0,-1),2,6,byrow = T)
rownames(con)<-c("Diet_3 vs Diet_5","Diet_3 vs Diet_6")
library(gmodels)
fit.contrast(model,"Treatment",con)

library(agricolae)
LSD.test(model,"Treatment") $group

#

con1<-c(0,0,1,0,-1,0)
library(gmodels)
fit.contrast(model,"Treatment",con1)

library(gmodels)
con2<-c(0,0,1,0,0,-1)
fit.contrast(model,"Treatment",con2)
#d
library(agricolae)
LSD.test(model,"Treatment")$group

##
#a significant difference between the treatments
blck<-factor(rep(paste0("Block_",1:5),each=6))
treat<-factor(rep(LETTERS[1:6],times=5))
b1<-c(33.5,38.9,40.9,41.9,41.9,42.8)
b2<-c(33,37.8,39.6,40.4,40.5,41.9)
b3<-c(32,37,39,39.9,39.8,40.7)
b4<-c(30,36,38.3,39.5,39.2,39.8)
b5<-c(28.5,34.8,37.5,38.7,38.8,39.4)
yld<-c(b1,b2,b3,b4,b5)
data<-data.frame(Block=blck,Treatment=treat,Yield=yld)
model<-aov(Yield~Block+Treatment, data=data)
summary(model)
#
con<-c(0,0,0,1,-1,0)
library(gmodels)
fit.contrast(model,"Treatment",con)
#95%
con<-c(0,1,-1,0,0,0)
library(gmodels)
fit.contrast(model,"Treatment",con, conf.int=0.95)
#
library(agricolae)
duncan.test(model,"Treatment")$group
LSD.test(model,"Treatment")$group

#row,c,t
 #a
r1<-c(192,190,214,221)
r2<-c(195,203,139,152)
r3<-c(292,218,245,204)
r4<-c(249,210,163,134)
y<-c(r1,r2,r3,r4)
row<-factor(rep(paste0("Cow-",1:4),each=4))
col<-factor(rep(paste0("Period-",1:4),4))
treat<-factor(paste0("Treat-",c(1:4,2,4,1,3,3,1,4,2,4,3,2,1)))
data<-data.frame(Row=row,Column=col,Treatment=treat,y=y)
model<-aov(y~Row+Column+Treatment, data=data)
summary(model)
#b
library(agricolae)
LSD.test(model,"Treatment")$group

#hit 
set.seed(125)
hmi<-function(a,b,g,n){
  x<-runif(n,a,b)
  k<-g(x)
  c<-max(k)
  y<-runif(n,0,c)
  v<-NULL
  for(i in 1:n){
    if(y[i]>k[i]){
      v[i]<-0
    }
    else{
      v[i]<-1
    }
  }
  nh<-sum(v)
  i<-(c*(b-a)*nh)/n
  cat("Hit or Miss Integration =",i)
}
g<-function(x){
  (sin(pi*x))^2
}
hmi(0,1,g,10000000)

#LCM
rand<-function(seed,a,c,m,n){
  vec<-NULL
  for( i in 1:n){
    z<-(a*seed+c)%%m
    seed<-z
    vec[i]<-z
  }
  data.frame ("Random number"=vec,"Between 0 and 1"=vec/m)
}
rand(52,21,53,100,5)
#CLCM
ran<-function(m,a,seed,n){
  k<-length(m)
  start<-rep(seed,k)
  vec<-c()
  for (i in 1:n){
    z<-rep(0,k)
    sum<-0
    for(j in 1:k){
      z[j]<-(a[j]*start[j])%%m[j]
      start[j]<-z[j]
      sum<-sum+((-1)^(j-1))*z[j]
    }
    x<-sum%%(m[1]-1)
    if(x>0){
      R<-x/m[1]
    } else{(m[1]-1)/m[2]}
    vec<-c(vec,R)
  }
  return(vec)
}
m<-c(2147483563,2147483399)
a<-c(40014,40692)
seed<-73
n<-10
ran(m,a,seed,n)

#control v
x<-function(f,g,rep,lo,up,control=T){
  mu<-integrate(f,lo,up)$value
  u<-runif(rep,lo,up)
  A<-g(u)
  B<-f(u)
  corr<-cor(A,B)
  a<--cov(A,B)/var(B)
  if(!control){
    T1<-A
  }else{
    T1<-A+a*(B-mu)
  }
  M<-(up-lo)*mean(T1)
  if(!control){
    list(MC=M)
  }else{
    list(correlation=corr,MC=M)
  }
}
f<-function(x){exp(-.5)/(1+x^2)}
g<-function(x){exp(-x)/(1+x^2)}
x(f,g,100000,0,1,control = T)



#Control Variates and Regression
x<-function(f,g,rep,lo,up){
  mu<-integrate(f,lo,up)$value
  u<-runif(rep,lo,up)
  L<-lm(g(u)~f(u))
  a<-L$coeff[2]
  theta.hat<-sum(L$coeff*c(1,mu))
  c(theta.hat,summary(L)$sigma^2,summary(L)$r.squared)               
}
f<-function(x){exp(-.5)/(1+x^2)}
g<-function(x){exp(-x)/(1+x^2)}
x(f,g,10000000,0,1)

#exp
rexpon<-function(size,lambda){ 
  X<-rep(0,size) 
  for(i in 1:size){
    U<-runif(1)
    X[i]<--(1/lambda)*(log(U))
  }
  X
}
set.seed(100)
x<- rexpon(10000,2)
x
mean(x)
var(x)
sd(x)

plot
library(MASS)
win.graph(6,3,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
truehist(x)
z<- (mean(x)-2)/(sd(x)/sqrt(10000))
#gamma
rgam<-function(size, shape, rate){ 
  x<-rep(0,size) 
  for ( i in 1:size){
    x[i]<-(-1/rate)*sum(log(runif(shape)))
  }
  x
}
set.seed(100)
x<- rgam(10000,3,2)
library(MASS)
win.graph(6,3,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
truehist(x)
 test the hypothesis that
z<- (mean(x)-(3/2))/sd(x)/(sqrt(10000))
#Box M
rnormal<-function(size, Mean, Sd){
  L<-size/2
  x<-rep(0,size)
  for(i in 1:L){
    u1<-runif(1)
    u2<-runif(1)
    x[2*i-1]<-Mean+Sd*sqrt(-2*log(u1))*cos(2*pi*u2)
    x[2*i]<-Mean+Sd*sqrt(-2*log(u1))*sin(2*pi*u2)
  }
  x
}
set.seed(100)
x<- rnormal(10000,5,2)
library(MASS)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
truehist(x)
# Ho:mu=5
mean(x)
sd(x)
z<- (mean(x)-5)/(sd(x)/sqrt(10000))
qnorm(0.005)
#Ho:sigma=4
sigma<-4
n<-10000
chi<-(n-1)*var(x)/sigma
qchisq(0.01,n-1)
#chi
rchisquare<-function(size, df){
  x<-rep(0,size)
  for ( i in 1:size){
    x[i]<-(-2)*sum(log(runif(df/2)))
  }
  x
}
set.seed(100)
x<- rchisquare (10000,10)
library(MASS)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
truehist(x)
lines (density(x))

#binom
rbino<-function(size,n,p){ 
  X<-rep(0, size)
  for(i in 1:size){
    x<-0
    sum<-choose(n,x)*(p^x)*(1-p)^(n-x)
    u<-runif(1)
    while(u>sum){
      x<-x+1
      sum<-sum+choose(n,x)*(p^x)*(1-p)^(n-x)
    }
    X[i]<-x
  }
  X
}
set.seed(100)
n<-7
p<-0.6
N<-10000
x<- rbino(N,n,p)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
freq<-table(x)
barplot(freq)
plot(freq)
mu<-n*p
phat<-mean(x)/n
z<-(mean(x)-mu)/sqrt(n*phat*(1-phat)/N)

#or
phat<- mean(x)/n
z<-(mean(x)/n-p)/sqrt(phat*(1-phat)/(n*N))

#rgeom

rgeometric<-function(size,prob){
  x<-rep(0,size)
  for(i in 1:size){
    u<-runif(1,0,1)
    x[i]<-floor(log(u)/log(1-prob))
  }
  x
}
set.seed(100)
x<- rgeometric (10000,0.25)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
freq<-table(x)
barplot(freq)
plot(freq,type="h")

#pois
rPoisson<-function(size,lambda){ 
  y<-rep(0, size) 
  for(i in 1:size){
    x<-0
    sum<-(exp(-lambda)*(lambda)^x)/factorial(x)
    u<-runif(1)
    while(u>sum){
      x<-x+1
      sum<-sum+(exp(-lambda)*(lambda)^x)/factorial(x)
    }
    y[i]<-x
  }
  y
}
set.seed(100)
N<-10000
lambda<-2
x<- rPoisson (N, lambda)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
freq<-table(x)
barplot(freq)
plot(freq,type="h")
lambhat<-mean(x)
z<-( lambhat-lambda)/sqrt(lambhat/N)

#marsaglia
rmarsag<-function(size, Mean, Sd){
  x<-rep(0,size)
  L<-size/2
  for(i in 1:L){
    n<-1
    while(n==1){
      u1<-runif(1)
      u2<-runif(1)
      w1<-(2*u1)-1
      w2<-(2*u2)-1
      w<-w1^2+w2^2
      if(w<1){
        x[2*i-1]<-Mean+Sd*sqrt(-2*log(w)/w)*(w1)
        x[2*i]<-Mean+Sd*sqrt(-2*log(w)/w)*(w2)
        n<-0
      }
    }
  }
  x
}
set.seed(100)
x<- rmarsag (10000,3,4)
library(MASS)
win.graph(6,2.5,8)
par(mfrow=c(1,2))
par(mar=c(4,4,1,1))
plot(density(x))
truehist(x)
lines (density(x))

#Gibs
set.seed(3)
rgibbs<-function(x0, y0, iter){
  x<-matrix(x0, iter)
  y<-matrix(y0, iter)
  for( i in 2:iter){
    u<-runif(1)
    x[i]<-sqrt(u*(1-x[i])^2)
    u<-runif(1)
    y[i]<-sqrt(u*(1-x[i])^2)
  }
  list(x=x, y=y)
}
r<-rgibbs(0.4,0.4,2000)
M<-cbind(x<-r$x,y<-r$y)
apply(M,2,mean)
cout("the estimates agree with their corresponding expected values")
#simulation of a poisson process
set.seed(3)
t<-0
ET<-NULL 
NF<-NULL
while(t<=2){ 
  u<-runif(1) 
  t1<-t-log(u)/5 
  if(t1>2){break} 
  ET<-c(ET,t1)
  NF<-c(NF,19+ceiling(21*u)) 
  t<-t+t1 
}
data.frame("Time"=ET, "Number of Fans"=NF)

From 0 to 0.3567 hour, 23 fans arrived (1st bus).
From 0.3567 to 0.3995 hour, 36 fans arrived (2nd bus)
##### 


lambda <- 3.8

# (a)
dpois(4, lambda)          # P(X = 4)
ppois(3, lambda)          # P(X < 4)
ppois(4, lambda)          # P(X ≤ 4)
1 - ppois(4, lambda)      # P(X > 4)
1 - ppois(3, lambda)      # P(X ≥ 4)

# (b)
ppois(7, lambda) - ppois(3, lambda)   # P(3 < X < 8)
ppois(8, lambda) - ppois(3, lambda)   # P(4 ≤ X ≤ 8)
ppois(7, lambda) - ppois(4, lambda)   # P(4 ≤ X < 8)
ppois(8, lambda) - ppois(4, lambda)   # P(4 < X ≤ 8)

# (c)
dpois(3, lambda) + dpois(4, lambda)   # P(X = 3 or X = 4)

# (d)
qpois(0.5, lambda)                    # Find 'a' such that P(X ≤ a) ≥ 0.5

# (e)
set.seed(123)
y <- rpois(70, lambda)                # Random sample of size 70
y
