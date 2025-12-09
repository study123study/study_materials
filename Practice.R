##############1
#T2 test
y1<-c(15,17,15, 13,20,15,15,13,14,17,15,17,15,18,18,15,18,10,18,18,13,16,11,16,16,18)
y2<-c(24,32 ,29,10, 26,26,26 ,22,30,30,26,28,29,32,31 ,26,33,19,30,34,30,16,25,26,23,34)
y3<-c(14,26,23,16,28,21 ,22,22,17 ,27 ,20,24,24,28,27,21 ,26,17 , 29,26,24,16,23,16,21,24)

hot<-matrix(c(y1,y2,y3),ncol=3)
hot
n<-nrow(hot)
n
p<-ncol(hot)
p
ybar<-colMeans(hot)
ybar
S<-var(hot)
S
mu0<-c(14,25,20)
T2<-n*t(ybar-mu0)%*%solve(S)%*%(ybar-mu0)
T2
ndf<-p
ddf<-n-p
fval<-qf(p=0.05,ndf,ddf,lower.tail = FALSE)
T2tab<-((n-1)*p/(n-p))*fval
 
#CI 
me<-sqrt(T2tab*diag(S/n))
me
cis<-cbind(ybar-me,ybar+me)
cis
#mean diff
med<-sqrt(T2tab*((S[1,1]+S[3,3]-2*S[1,3])/n))
med
cid<-cbind(ybar[1]-ybar[3]-med,ybar[1]-ybar[3]+med)
cid

########2

n1<-271
n2<-138
n3<-107
n<-n1+n2+n3
n
p<-4
g<-3
x1b<-matrix(c(2.066,0.480,0.082,0.360),ncol=1)
x1b
x2b<-matrix(c(2.167,0.596,0.124,0.418),ncol=1)
x2b
x3b<-matrix(c(2.273,0.521,0.125,0.383),ncol=1)
x3b
s1<-matrix(c(0.291,-0.001,0.002,0.010,-0.001,0.011,0.000,0.003,0.002,0.000,0.001,0.000,0.010,0.003,0.000,0.010),ncol=4)
s1
s2<-matrix(c(0.561,0.011,0.001,0.037,0.011,0.025,0.004,0.007,0.001,0.004,0.005,0.002,0.037,0.007,0.002,0.019),ncol=4)
s2
s3<-matrix(c(0.261,0.030,0.003,0.018,0.030,0.017,-0.000,0.006,0.003,-0.000,0.004,0.001,0.018,0.006,0.001,0.013),ncol=4)
s3

w<-(n1-1)*s1+(n2-1)*s2+(n3-1)*s3
w
xb<-((n1*x1b)+(n2*x2b)+(n3*x3b))/n
xb
b<-(n1)*(x1b-xb)%*%t(x1b-xb)+(n2)*(x2b-xb)%*%t(x2b-xb)+(n3)*(x3b-xb)%*%t(x3b-xb)
b
wilk<-det(w)/det(b+w)
wilk

# for p>=1 and g=3
fcal<-((n-p-2)/p)*((1-sqrt(wilk))/sqrt(wilk))
ftab<-qf(0.01,2*p,2*(n-p-2),lower.tail = FALSE)
cat("Calculated value =",fcal,"\n")
cat("Tabulated vaue =",ftab, "\n")

# for large n=516 the exact test
chcal<--(n-1-(p+g)/2)*log(wilk)
chtab<-qchisq(0.01,p*(g-1),lower.tail = FALSE)
cat("Calculated value =",chcal,"\n")
cat("Tabulated vaue =",chtab, "\n")

t1h<-x1b-xb
t2h<-x2b-xb
t3h<-x3b-xb
a<-t1h[3,1]-t3h[3,1]
b<-t1h[3,1]-t2h[3,1]
c<-t2h[3,1]-t3h[3,1]
crv<-qt((0.05/(p*g*(g-1))),n-g,lower.tail = FALSE)
d<-sqrt((w[3,3]/(n-g))*((1/n1)+(1/n3)))
e<-sqrt((w[3,3]/(n-g))*((1/n1)+(1/n2)))
f<-sqrt((w[3,3]/(n-g))*((1/n2)+(1/n3)))

cat("CI for tau13-tau33 =",a-crv*d,a+crv*d,"\n","\n")
cat("CI for tau13-tau23 =",b-crv*e,b+crv*e,"\n","\n")
cat("CI for tau23-tau33 =",c-crv*f,c+crv*f,"\n","\n")

sp<-w/((n1-1)+(n2-1)+(n3-1))
m<-((n1-1)+(n2-1)+(n3-1))*log(det(sp))-((n1-1)*log(det(s1))+(n2-1)*log(det(s2))+(n3-1)*log(det(s3)))
m
u<-(1/(n1-1)+1/(n2-1)+1/(n3-1)-1/((n1-1)+(n2-1)+(n3-1)))*((2*p*p+3*p-1)/(6*(p+1)*(g-1)))
u
chcal<-(1-u)*m
chtab<-qchisq(0.05,(p*(p+1)*(g-1)/2),lower.tail = FALSE)
cat("Calculated value =",chcal,"\n","\n")
cat("Tabulated vaue =",chtab, "\n")

#3
skull<-read.table("F:/Statistics semester wise book & sheet/7th semester/STA 4106/RCODE/T6-13.dat")

colnames(skull)<-c("MaxBreath","BasHeight","BasLength","NasHeight","TunePeriod")

head(skull)

x<-as.matrix(skull[,1:4])
time<-skull[,5]

n1<-sum(time==1)
n2<-sum(time==2)
n3<-sum(time==3)
n<-n1+n2+n3
ens<-c(n1,n2,n3)
p<-4
g<-3

s1<-cov(x[time==1,])
s2<-cov(x[time==2,])
s3<-cov(x[time==3,])

W<-(n1-1)*s1+(n2-1)*s2+(n3-1)*s3
spooled<-W/(n-g)

# Test equality of the Covariance Matrices or Box's M test
M<-(n-g)*log(det(spooled))-((n1-1)*log(det(s1))+(n2-1)*log(det(s2))+(n3-1)*log(det(s3)))
u<-(sum(1/(ens-1))-1/sum(ens-1))*(2*p*p+3*p-1)/(6*(p+1)*(g-1))
c<-(1-u)*m
df<-p*(p+1)*(g-1)/2
pval<-pchisq(c,df,lower.tail = FALSE)                                          # 0.3942866

xb=colMeans(x)
t1h<-colMeans(x[time==1,])-xb
t2h<-colMeans(x[time==2,])-xb
t3h<-colMeans(x[time==3,])-xb
th<-cbind(t1h,t2h,t3h)

B<-th%*%diag(ens)%*%t(th)
wilks<-det(W)/det(W+B)                                                         # 0.8301027
bart<--(n-1-((p+g)/2))*log(wilks)
asy.pval<-pchisq(bart,p*(g-1),lower.tail = FALSE)                              # 0.04353074


# 'matplot' plot all the columns of one matrix against all the columns of another.
# Here we plot the 4 columns of tauhat against the one column c(1,2,3,4).

matplot(1:4, th, type='l',xaxp=c(1,4,3),lty=c(1,2,4),col=c(1,2,1),xlab="Variables",ylab="Treatment Effects")
legend("topright",c("tau1","tau2","tau3"),lty=c(1,2,4),col=c(1,2,1))

fit<-manova(x~as.factor(time))   # Note: time is a vector of labels - 'factor'
summary(fit)
summary(fit, test = "Wilks")

#Critical value F with df 2p and 2*(n-p-2)
#This is the Distribution of wilks' Lambda because p=4, g=3 
crit.value <- qf (0.95, 8, 168)
crit.value

# Conclusion: with F test statistic = 2.049 > 1.9939, and a p-value of 0.04358<0.05 = ??,
# we reject our HO and conclude that the time effect differences exist.
# There is a difference of male Egyptian skulls for three different time periods.

summary(fit, test = "Pillai")
summary(fit, test = "Hotelling-Lawley")
summary(fit, test = "Roy") 
summary.aov(fit)

#Bonferroni intervals on all differences

alpha = 0.05    # Some are significant at alpha = 0.05 m = p*g*(g-1)/2
m<-p*g*(g-1)/2
q = qt(alpha/(2*m), n-g, lower.tail = FALSE)
confint<-array(dim=c(g,g,p,3),dimnames = list(NULL, NULL, NULL, c("lower", "point", "upper")))

for (k in 1:g) {
  for (l in 1:g) {
    halfwidth = q*sqrt((diag (w)/(n-g))*(1/ens [k] + 1/ens [1]))
    mid = th [, k] - th[,1]
    lower = mid - halfwidth
    upper = mid + halfwidth
    confint [k, 1,, ] = cbind(lower, mid, upper)
  }
}
for(k in 2:g){
  for (l in 1:(k-1)) {
    mat=round(confint[k,l,,],3)
    sig=sign(mat[,1]*mat[,3])
    sig[sig==-1]=''
    sig[sig==1]='*'
    mat=cbind(mat,sig)
    cat(paste(100*(1-alpha),"% Bonferrani CI on tau",k,"minus tau",l,"are:","\n"))
    prmatrix(mat,mat,quote = FALSE)
    cat("\n")
  }
}

# pair comparison
g<-3
p<-4

#All Tune Periods are 30
n1<-length (which((skull$TunePeriod ==1)))
n2<-length (which((skull$TunePeriod ==2)))
n3<-length (which((skull$TunePeriod ==3)))
n<-n1+n2+n3

#Finding means across the columns and we don't want Tune Period so we remove it!
x1b<-colMeans(skull[skull$TunePeriod==1, -5])
x2b<-colMeans(skull[skull$TunePeriod==2,-5])
x3b<-colMeans(skull[skull$TunePeriod==3, -5])
xb<-(n1*x1b+n2*x2b+n3*x3b)/(n1+n2+n3)

#Finding SE
S1<-cov(skull[skull$TunePeriod==1, -5])
S2<-cov(skull[skull$TunePeriod==2, -5])
S3<-cov(skull[skull$TunePeriod==3, -5])
W<-(n1-1)*S1+(n2-1)*S2+(n3-1)*S3

alpha<-0.05
crit.b<-qt(1-alpha/(p*g*(g-1)),n-g)

for(i in 1:p){
  LCI12<-(x1b[i]-x2b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n2))
  UCI12<-(x1b[i]-x2b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n2))
  cat("mu1[",i,"]-mu2[",i,"]=(",LCI12,",",UCI12,")\n",sep="")
  
  #\mu_{11}-\mu_{31}
  LCI13<-(x1b[i]-x3b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n3))
  UCI13<-(x1b[i]-x3b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n1+1/n3))
  cat("mu1[",i,"]-mu3[",i,"]=(",LCI13,",",UCI13,")\n",sep="")
  
  #\mu_{21}-\mu_{31}
  LCI23<-(x2b[i]-x3b[i])-crit.b*sqrt(W[i,i]/(n-g)*(1/n2+1/n3))
  UCI23<-(x2b[i]-x3b[i])+crit.b*sqrt(W[i,i]/(n-g)*(1/n2+1/n3))
  cat("mu2[",i,"]-mu3[",i,"]=(",LCI23,",",UCI23,")\n",sep="")
}

# Test Normality
#Time Period 1
S1inv<-solve(S1) # Inverse the matrix s1. 
skull1<- skull[skull$TunePeriod==1, -5]
datachisq<-diag(t(t (skull1)-x1b)%*%S1inv%*%(t(skull1)-x1b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")

#Time Period 2
S2inv<-solve(S2) # Inverse the matrix s2. 
skull2<- skull[skull$TunePeriod==2, -5]
datachisq<-diag(t(t (skull2)-x2b)%*%S2inv%*%(t(skull2)-x2b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")

#Time Period 3
S3inv<-solve(S3) # Inverse the matrix s3. 
skull3<- skull[skull$TunePeriod==3, -5]
datachisq<-diag(t(t (skull3)-x3b)%*%S3inv%*%(t(skull3)-x3b))
qqplot(qchisq(ppoints(500),df=p),datachisq, main="",xlab="Theoretical Quantiles", ylab="Sample Quantiles")


#######4


z0<-rep(1,5)
z1<-c(0,1,2,3,4)
y1<-c(1,4,3,8,9)
y2<-c(-1,-1,2,3,2)

y<-matrix(c(y1,y2),ncol=2)
y

z<-matrix(c(z0,z1),ncol=2)
z
tz.z<-t(z)%*%z
tz.z

tz.z.inv<-solve(tz.z)
tz.z.inv

tz.y1<-t(z)%*%y1
tz.y1

b1<-solve(t(z)%*%z)%*%t(z)%*%y1
b1

tz.y2<-t(z)%*%y2
tz.y2

b2<-solve(t(z)%*%z)%*%t(z)%*%y2
b2

b<-matrix(c(b1,b2),ncol=2)
b

y.hat<-z%*%b
y.hat

e.hat<-(y-y.hat)
e.hat


# Properties Check
sum(e.hat)
t(y.hat)%*%e.hat
t(z)%*%e.hat

ssr<-t(y.hat)%*%y.hat
ssr

sse<-t(e.hat)%*%e.hat
sse

sst<-t(y)%*%y
sst
ssr+sse


# By Regreesion model
m<-lm(y~z1)
summary(m)


