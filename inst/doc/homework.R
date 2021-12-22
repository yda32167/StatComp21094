## -----------------------------------------------------------------------------
x=3
sqrt(x)

## -----------------------------------------------------------------------------
plot(c(1:6),c(1:6),main="test",xlim = c(0,7),ylim = c(0,7))

## -----------------------------------------------------------------------------
cityofstudent= c("武汉","合肥","南京","合肥","上海")
table(cityofstudent)

## -----------------------------------------------------------------------------
set.seed(123)
n<-1000
u<-runif(n)
y<-seq(0,20,.01)

par(mfrow=c(1,2))
#sigma=2的情形，可以看出拟合较好
sig_1<-2


x_1<-sqrt(-2*sig_1^2*log(1-u))
hist(x_1,probability = T,main=expression(sigma==2))
#添加理论密度曲线,图中为红线
lines(y,y/sig_1^2*exp(-y^2/(2*sig_1^2)),col="red")

#sigma=5的情形，可以看出拟合较好
sig_2<-5

x_2<-sqrt(-2*sig_2^2*log(1-u))
hist(x_2,probability = T,main=expression(sigma==5))
#添加理论密度曲线,图中为红线

lines(y,y/sig_2^2*exp(-y^2/(2*sig_2^2)),col="red")


## -----------------------------------------------------------------------------
bif<-function(p)
{
#定义函数bif,输入p_1=p,输出对应混合分布直方图
n<-1e3
X1<-rnorm(n,0,1)
X2<-rnorm(n,3,1)
p_1<-p;p_2<-1-p_1
r<-sample(c(0,1),n,replace = T,prob = c(p_2,p_1))
Z<-r*X1+(1-r)*X2
hist(Z,main="")
}

## -----------------------------------------------------------------------------
set.seed(123)
#输出p_1=0.75时的混合分布直方图
bif(0.75)
title(main="p=0.75",outer = T)
#输出p_1=0.1,...,1时的混合分布直方图
par(mfrow=c(2,5))

## -----------------------------------------------------------------------------
n<-1e4
lambda1<-5
t<-10
Nt1<-rpois(n,lambda1*t)
afa<-1
beta<-1
#此处将Yi的概率密度参数设置为Ga（1，1）
Xt1<-rgamma(n,afa*Nt1,beta)
#由于Yi独立同分布，其和服从Ga（Nt1，1）
hist(Xt1,prob=TRUE)
EX1<-sum(Xt1)*{1/n}
VarX1<-sum((Xt1-EX1)^2)*{1/(n-1)}
#此处计算出由随机数生成的样本观测值的均值和方差
TRUEEX1<-t*lambda1*afa/beta
TRUEVarX1<-t*lambda1*(afa^2+afa)/beta
#此处计算出由公式得出的理论均值和方差
EX1/TRUEEX1
VarX1/TRUEVarX1
#通过相除得出λ=5时，理论与实际比值在1附近
lambda2<-10
Nt2<-rpois(n,lambda2*t)
Xt2<-rgamma(n,afa*Nt2,beta)
hist(Xt2,prob=TRUE)
EX2<-sum(Xt2)*{1/n}
VarX2<-sum((Xt2-EX2)^2)*{1/(n-1)}
TRUEEX2<-t*lambda2*afa/beta
TRUEVarX2<-t*lambda2*(afa^2+afa)/beta
EX2/TRUEEX2
VarX2/TRUEVarX2
#此处取λ=10，比值依旧在1附近

## -----------------------------------------------------------------------------
g<-function(x){mean((1/beta(3,3))*(x^2)*(1-x)^2)}
m <- 10000#重复10000次
estimates=matrix(0,9,2)
for (i in 1:9) {
  x <- runif(m,0,i/10)#将9次试验的值赋值到矩阵里
    estimates[i, 1] <- (i/10)*g(x)
    estimates[i, 2] <- pbeta(i/10,3,3)}
estimates#展示函数的估计值与真实值，第一列为估计值，第二列为真实值

## -----------------------------------------------------------------------------
MC.Phi <- function(x,s, R = 10000, antithetic = TRUE) {
u <- runif(R/2)
if (!antithetic) v <- runif(R/2) else
v <- 1 - u
u <- c(u, v)
cdf <- numeric(length(x))
for (i in 1:length(x)) {se
g <- 1-exp(-(u * x[i])^2 / (2*s^2))
cdf[i] <- mean(g)
}
cdf}#定义对偶函数
m <- 1000
MC1 <- MC2 <- numeric(m)
x <- 2
se=matrix(0,10,1)#对方差进行不同取值，从1，2到10依次选取
for (j in 1:10) {
for (i in 1:m) {
MC1[i] <- MC.Phi(x,j, R = 1000, anti = FALSE)
MC2[i] <- MC.Phi(x,j, R = 1000)
}
se[j,1]=((var(MC1) - var(MC2))/var(MC1))
}
se

## -----------------------------------------------------------------------------
m <- 100000
theta.hat <- se <- numeric(2)
g <- function(x) {(x^2)/(sqrt(2*pi))*exp(-x^2/2) * (x > 1)}
x <- rexp(m, 1)+1 #第一个函数取近似的指数函数，取值范围为x>1
fg <- g(x) / exp(-x+1)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
x <- rnorm(m) #第二个函数取标准正态函数
fg <- g(x) / dnorm(x)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
rbind(theta.hat, se)

## -----------------------------------------------------------------------------
m <- 100000
g <- function(x) {(x^2)/(sqrt(2*pi))*exp(-x^2/2) * (x > 1)}
x <- rexp(m, 1)+1 #函数取近似的指数函数，取值范围为x>1
fg <- g(x) / exp(-x+1)
mean(fg)


## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 2
m <- 1000 
uci<-lci <- numeric(m)
c=0
for (i in 1:m) {
x <- rchisq(n,df=2)
uci[i] <- mean(x)+qt(1-alpha/2,n-1)*sd(x)/sqrt(n-1)
lci[i] <- mean(x)-qt(1-alpha/2,n-1)*sd(x)/sqrt(n-1)
c=c+sum(lci[i]<2&uci[i]>2)
}
print(c/10000)

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rnorm(n, mean = 0, sd = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
M<-mean(UCL > 4)
M

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- rchisq(n,df=1)
ttest <- t.test(x, alternative = "two.sided", mu = mu0)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- runif(n,0,2)
ttest <- t.test(x, alternative = "two.sided", mu = mu0)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
mu0 <- 1
m <- 10000 #number of replicates
p <- numeric(m) #storage for p-values
for (j in 1:m) {
x <- rexp(n,1)
ttest <- t.test(x, alternative = "two.sided", mu = mu0)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
se.hat <- sqrt(p.hat * (1 - p.hat) / m)
print(c(p.hat, se.hat))

## -----------------------------------------------------------------------------
library(MASS)
n <- c(10, 20,50,100) #样本大小
d <- 4#多元正态的维数
cv1 <- qchisq(.95, d*(d+1)*(d+2)/6) #确定维数下的拒绝域上界临界值（还需乘以6/n）

sk <- function(x) {
#计算样本偏度系数
mu<-numeric(d)
sigma<-diag(d)
x <- mvrnorm(n[i],mu,sigma)
xbar<-colMeans(x)
SIGMA<- matrix(data=0, nrow = d, ncol =d)
for (a in 1:n[i]){
  SIGMA=SIGMA+(1/n[i])*(x[a,]-xbar)%o%(x[a,]-xbar)
}
SIGMA1=solve(SIGMA)
beta1=0
for (a in 1:n[i]){
  for (b in 1:n[i]){
    beta1=((1/n[i])^2)*((x[a,]-xbar)%*%SIGMA1%*%(x[b,]-xbar))^3+beta1
}
}
return( beta1)
}
p.reject <- numeric(length(n)) #保存样本偏度系数的模拟结果
m <- 200 #每次实验重复200次
for (i in 1:length(n)) {
sktests <- numeric(m)
for (j in 1:m) {
#拒绝原假设为1
sktests[j] <- as.integer(abs(sk(x)) >= 6*cv1/n[i])
}
p.reject[i] <- mean(sktests) #拒绝的概率
}
p.reject

## -----------------------------------------------------------------------------
alpha <- .1
d=5
n <- 30
m <- 100
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N)
sk <- function(x) {
#计算样本偏度系数
mu<-numeric(d)
sigma <- diag(c(sample(c(1, 10), replace = TRUE,
size =d, prob = c(1-e, e))))
x <- mvrnorm(n,mu,sigma)
xbar<-colMeans(x)
SIGMA<- matrix(data=0, nrow = d, ncol =d)
for (a in 1:n){
  SIGMA=SIGMA+(1/n)*(x[a,]-xbar)%o%(x[a,]-xbar)
}
SIGMA1=solve(SIGMA)
beta1=0
for (a in 1:n){
  for (b in 1:n){
    beta1=((1/n)^2)*((x[a,]-xbar)%*%SIGMA1%*%(x[b,]-xbar))^3+beta1
}
}
return( beta1)
}
cv <- qchisq(.95, d*(d+1)*(d+2)/6) #确定维数下的拒绝域上界临界值（还需乘以6/n）
for (j in 1:N) { 
e <- epsilon[j]
sktests <- numeric(m)
for (i in 1:m) { 
sktests[i] <- as.integer(abs(sk(x)) >= 6*cv/n)
}
pwr[j] <- mean(sktests)
}
#画出图形
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) 
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
library(MASS)
x<-scor
la<-function(x){
  xbar<-colMeans(x)
  n=88
SIGMA<- matrix(data=0, nrow = 5, ncol =5)
for (a in 1:n){
  t<-as.numeric(x[a,]-xbar)
  SIGMA=SIGMA+(1/n)*(t)%o%(t)
}
ev<-eigen(SIGMA)
lambda<-sum(ev$val)
lambda1<-ev$val[1]
theta<-lambda1/lambda
return(theta)
}
lambda.hat<-la(x)

B <- 1e2; 
thetastar <- numeric(B)
for(b in 1:B){
mb <- 1:88
nmb <- sample(mb,replace=TRUE)
xstar=matrix(data=0, nrow = 88, ncol =5)
for (j in 1:88) {
xstar[j,]<-as.numeric(x[nmb[j],])
}
thetastar[b] <-la(xstar) 
}

round(c(bias=mean(thetastar)-lambda.hat,
se=sd(thetastar)),2)

## -----------------------------------------------------------------------------
library(bootstrap)
library(MASS)
x<-scor
la<-function(x){
  xbar<-colMeans(x)
  n=dim(x)[1]
SIGMA<- matrix(data=0, nrow = 5, ncol =5)
for (a in 1:n){
  t<-as.numeric(x[a,]-xbar)
  SIGMA=SIGMA+(1/n)*(t)%o%(t)
}
ev<-eigen(SIGMA)
lambda<-sum(ev$val)
lambda1<-ev$val[1]
theta<-lambda1/lambda
return(theta)
}
lambda.hat<-la(x)
thetastar<-numeric(88)
for(b in 1:88){
xt<-x[-b,]
thetastar[b] <-la(xt) 
}

round(c(bias=87*(mean(thetastar)-lambda.hat),
se=sd(thetastar)),2)

## -----------------------------------------------------------------------------
library(bootstrap)
library(MASS)
x<-scor
la<-function(x){
  xbar<-colMeans(x)
  n=88
SIGMA<- matrix(data=0, nrow = 5, ncol =5)
for (a in 1:n){
  t<-as.numeric(x[a,]-xbar)
  SIGMA=SIGMA+(1/n)*(t)%o%(t)
}
ev<-eigen(SIGMA)
lambda<-sum(ev$val)
lambda1<-ev$val[1]
theta<-lambda1/lambda
return(theta)
}
lambda.hat<-la(x)

B <- 1e2; 
thetastar <- numeric(B)
for(b in 1:B){
mb <- 1:88
nmb <- sample(mb,replace=TRUE)
xstar=matrix(data=0, nrow = 88, ncol =5)
for (j in 1:88) {
xstar[j,]<-as.numeric(x[nmb[j],])
}
thetastar[b] <-la(xstar) 
}
pc=lambda.hat+qnorm(0.95)*sd(thetastar)
BCa=c(lambda.hat-qnorm(0.975)*sd(thetastar),lambda.hat+qnorm(0.975)*sd(thetastar))
l=BCa[1]
u=BCa[2]
round(c(percentile=pc,lowerBCa =l,upperBCa =u),3)

## -----------------------------------------------------------------------------
#standard normal
n=1000
cv1 <- qnorm(.975) 
cv2 <- qnorm(.025) 
p.reject<-X <- numeric(n) #保存样本偏度系数的模拟结果
for (i in 1:n) {
X[i]=rnorm(1,0)
}
xbar<-mean(X)
xse<-sd(X)
for (j in 1:n) {
#拒绝原假设为1
p.reject[j] <- as.integer(X[j]>=cv1)|as.integer(X[j]<=cv2 )
}
p.reject1=mean(p.reject)
#standard normal sampling confidence interval
sv1 <- xbar+qnorm(.975) *xse
sv2 <-xbar+ qnorm(.025)*xse 
for (j in 1:n) {
#拒绝原假设为1
p.reject[j] <- as.integer(X[j]>=sv1)|as.integer(X[j]<=sv2 )
}
p.reject2=mean(p.reject)
#standard normal bootstrap confidence interval
for (i in 1:n) {
X[i]=rnorm(1,0)
}
thetastar <- numeric(1000)
for(b in 1:1000){
xstar <- sample(X,replace=TRUE)
thetastar[b] <- mean(xstar)
}
bv1 <-mean(thetastar)+qnorm(.975)*sd(thetastar)
bv2 <-mean(thetastar)+qnorm(.025)*sd(thetastar)
for (j in 1:1000) {
#拒绝原假设为1
p.reject[j] <- as.integer(thetastar[j]>=bv1)|as.integer(thetastar[j]<=bv2 )
}
p.reject3=mean(p.reject)
round(c(p1=p.reject1,p2=p.reject2,p3=p.reject3),3)

## -----------------------------------------------------------------------------
#standard normal
n=1000
cv1 <- qchisq(.975,5) 
cv2 <- qchisq(.025,5) 
p.reject<-X <- numeric(n) #保存样本偏度系数的模拟结果
X=rchisq(df=5,n)
xbar<-mean(X)
xse<-sd(X)
for (j in 1:n) {
#拒绝原假设为1
p.reject[j] <- as.integer(X[j]>=cv1)|as.integer(X[j]<=cv2 )
}
p.reject1=mean(p.reject)
#standard normal sampling confidence interval
sv1 <- xbar+qchisq(.975,5) *xse
sv2 <-xbar+ qchisq(.025,5)*xse 
for (j in 1:n) {
#拒绝原假设为1
p.reject[j] <- as.integer(X[j]>=sv1)|as.integer(X[j]<=sv2 )
}
p.reject2=mean(p.reject)
#standard normal bootstrap confidence interval
X=rchisq(df=5,n)
thetastar <- numeric(1000)
for(b in 1:1000){
xstar <- sample(X,replace=TRUE)
thetastar[b] <- mean(xstar)
}
bv1 <-mean(thetastar)+qchisq(.975,5)*sd(thetastar)
bv2 <-mean(thetastar)+qchisq(.025,5)*sd(thetastar)
for (j in 1:1000) {
#拒绝原假设为1
p.reject[j] <- as.integer(thetastar[j]>=bv1)|as.integer(thetastar[j]<=bv2 )
}
p.reject3=mean(p.reject)
round(c(p1=p.reject1,p2=p.reject2,p3=p.reject3),3)

## -----------------------------------------------------------------------------
    set.seed(212)
    N <- 5000               #length of chain
    burn <- 1000            #burn-in length
    X <- matrix(0, N, 2)    #the chain, a bivariate sample
    n1 <- 200                 #parameter n
    a <- 2
    b <- 3

    ###### generate the chain #####

    X[1, ] <- c(0.5, 100)            #initialize

    for (i in 2:N) {
        x2 <- X[i-1, 2]
        X[i, 1] <- rbeta(1,x2+a,n1-x2+b)
        x1 <- X[i, 1]
        X[i, 2] <- rbinom(1,n1,x1)
    }

    b <- burn + 1
    x <- X[b:N, ]

    # compare sample statistics to parameters
    colMeans(x)
    cov(x)
    cor(x)
    jchs<-t(x)
    for (v in 1:n1){
    jchs[v] <-  jchs[v]+10 
    jchs[v+1] <-  jchs[v]+10
    jchs[2*v] <-  jchs[v]+20 
    if (jchs[2*v]>jchs[v])
        jchs[v]<-1
    else jchs[v]<-0
    }
    
    plot(x, main="", cex=.5, xlab=bquote(X[1]),
         ylab=bquote(X[2]), ylim=range(x[,2]))
    X<- t(X)
    Gelman.Rubin <- function(psi) {
        # psi[i,j] is the statistic psi(X[i,1:j])
        # for chain in i-th row of X
        psi <- as.matrix(psi)
        n <- ncol(psi)
        k <- nrow(psi)

        psi.means <- rowMeans(psi)     #row means
        B <- n * var(psi.means)        #between variance est.
        psi.w <- apply(psi, 1, "var")  #within variances
        W <- mean(psi.w)               #within est.
        v.hat <- W*(n-1)/n + B/n+(B/(n*k))     #upper variance est.
        r.hat <- v.hat / W             #G-R statistic
        return(r.hat)
        }

    #plot the sequence of R-hat statistics
    rhat <- rep(0, N)
    for (j in b:N){
        rhat[j] <- Gelman.Rubin(X[,1:j])}
    plot(rhat[(b+1):N], type="l", xlab="", ylab="R",ylim=c(1,10))
    abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
#(1)第k项的计算
cpe<-function(d,a,k){
  g=((-1)^k*gamma((d+1)/2))/((2*k+1)*(2*k+2))*exp((2*k+2)*log(sqrt(sum(a^2)))+lgamma(k+3/2)-lgamma(k+1)-k*log(2)-lgamma(k+d/2+1))#第k项的计算公式
g#返回值
  }

## -----------------------------------------------------------------------------
#(2)给定k后，和的计算
cpe<-function(a,k,e){
  sm=0
  g0=0
  d<-length(a)
  for (i in 0:k) {
    g=((-1)^i*gamma((d+1)/2))/((2*i+1)*(2*i+2))*exp((2*i+2)*log(sqrt(sum(a^2)))+lgamma(i+3/2)-lgamma(i+1)-i*log(2)-lgamma(i+d/2+1))#第k项的计算公式
    sm=g+sm
    if(abs(g0+g)<e)break#如果两次计算的和小于给定值，就跳出循环
    else go=g
  }
  N=1000
  jchs<-matrix(0, N, 2)
    for (v in 1:100){
    jchs[v] <-  jchs[v]+10 
    jchs[v+1] <-  jchs[v]+10
    jchs[2*v] <-  jchs[v]+20 
    if (jchs[2*v]>jchs[v])
        jchs[v]<-1
    else jchs[v]<-0
    }
  m=c(i,sm)#输出数据为迭代次数以及所求公式的和
  return(m)
}

#(3)给定a为(1,2)时，对于不同k，e的取值进行试验
#k=5，e=1e-10
cpe(a=c(1,2),5,1e-10)

#k=10,e=1e-10
cpe(a=c(1,2),10,1e-10)

#k=20,e=1e-10
cpe(a=c(1,2),20,1e-10)

#k=100,e=1e-100
cpe(a=c(1,2),100,1e-100)

## -----------------------------------------------------------------------------
#第十四题公式
slu4 <- function(a, k){
  v1 =sqrt(a^2*(k-1)/(k-a^2))
  v2 =sqrt(a^2*k/(k+1-a^2))
  pt(v1, df = k-1) - pt(v2, df = k)
}
#第一步，先画出两函数差的图像：左减右
f<-function(u,k){
  (1+u^2/k)^(-(k+1)/2)#定义被积函数
}
slu5=function(a,k){#定义两函数的差
  v1 =sqrt(a^2*(k-1)/(k-a^2))
  v2 =sqrt(a^2*k/(k+1-a^2))
  gg=(log(k) - log(k-1))/2 + (2*lgamma(k/2) - lgamma((k+1)/2)- lgamma((k-1)/2)) + log(integrate(f,lower = 0,upper =v1 ,rel.tol =1e-5,k=k-1)$value)-log(integrate(f,lower = 0,upper =v2 ,rel.tol =1e-5,k=k)$value)
  return(gg)
}
k=4
#k=4时的图像
a=seq(-2,2,0.01)#在如下区间里进行作图
v=rep(0,length(a))
for (i in 1:length(a)) {
  v[i]=slu5(a[i],k)
}
plot(a,v,type = "l")
#可以看到在接近1.5的地方，函数有根。

#（2）下面求根
it <- 0
eps <- 1e-50#误差阈值
vec=c(4:25,100,500,1000)
sluo5 <- sapply(vec, function(k){uniroot(slu5, interval = c(1, 2), k=k)$root})
sluo4 <- sapply(vec, function(k){uniroot(slu4, interval = c(1, 2), k=k)$root})
root <- cbind(sluo5, sluo4)
colnames(root) <- c("11.5的根", "11.4的根")
rownames(root) <- as.character(vec)

knitr::kable(root)

## -----------------------------------------------------------------------------
lambda_hat=1
for (i in 1:9) {
  lambda_hat<- (7*exp(lambda_hat)-7)/(exp(lambda_hat)+9.8)
}
lambda_hat

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
zt1<-lapply(formulas, lm,data=mtcars)
zt2<-lapply(zt1,rsq)
print(zt2)

## -----------------------------------------------------------------------------
x<-data.frame(cbind(x1=seq(2,8,2),x2=c(1:4)))
zt<-vapply(x, sd, FUN.VALUE = c('a'=0))
print(zt)
#结果显示第一列向量（2，4，6，8）的标准差为2.581989，第二列向量（1，2，3，4）的标准差为1.290994

## -----------------------------------------------------------------------------
x<-data.frame(a=1:5,b=runif(5),c=c(TRUE,FALSE,TRUE,FALSE,TRUE))
zt1<-vapply(x, is.numeric, logical(1))
zt2<-vapply(x[zt1], sd, FUN.VALUE = c('a'=0))
print(zt2)

## -----------------------------------------------------------------------------

function (cl = NULL, X, fun, ..., chunk.size = NULL) 
{
    cl <- defaultCluster(cl)
    nchunks <- staticNChunks(length(X), length(cl), chunk.size)
    do.call(c, clusterApply(cl = cl, x = splitList(X, nchunks), 
        fun = sapply, FUN = fun, ...), quote = TRUE)
}




## -----------------------------------------------------------------------------
library(png)
ans1 <- readPNG("~/StatComp21094/data/pic/tuone.png")
r <- nrow(ans1)/ncol(ans1) 
plot(c(0,1),c(0,r))
rasterImage(ans1,0,0,1,r)

