## -----------------------------------------------------------------------------
rm(list=ls())
library(SA24204187)
set.seed(123)
random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
CG.mean(random_matrix_normal)
selectCG.mean(random_matrix_normal,0.05,20)

## -----------------------------------------------------------------------------
rm(list=ls())
library(SA24204187)
set.seed(123)
random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
trimmed.mean(random_matrix_normal)

## -----------------------------------------------------------------------------
rm(list=ls())
library(SA24204187)
set.seed(123)
random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
robust.mean(random_matrix_normal,method = "huber") 

## -----------------------------------------------------------------------------
rm(list=ls())
library(SA24204187)
set.seed(123)
random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
MOM.mean(random_matrix_normal,k = 40,method = "geometric") 

## -----------------------------------------------------------------------------
rm(list=ls())
library(SA24204187)
set.seed(123)
random_matrix_normal <- matrix(rnorm(20000, mean = 0, sd = 1), nrow = 200, ncol = 100)
outliers = matrix(rnorm(2000, mean = 5, sd = 10), nrow = 200, ncol = 10)
X <- cbind(random_matrix_normal,outliers)
epsilon <- 0.1
# 调用投影梯度下降算法
w <- projected_gradient_descent(X, epsilon)
# 计算加权均值
mu_w <- X %*% w
# 输出结果
print("加权均值:")
print(mu_w)

## -----------------------------------------------------------------------------
attach(iris)
head(iris)

## -----------------------------------------------------------------------------
table(iris$Species) #Check the number of flowers in different species
library(xtable)
xtable::xtable(table(iris$Species))
detach(iris)

## -----------------------------------------------------------------------------
plot(freeny.x[,1],freeny.x[,2],xlab = "lag quarterly revenue",ylab = "price index")
LM = lm(freeny.x[,1]~freeny.x[,2])
summary(LM)

## -----------------------------------------------------------------------------
head(swiss)

## -----------------------------------------------------------------------------
plot(swiss)


## -----------------------------------------------------------------------------
rm(list=ls())
n <- 1000
u <- runif(n) # Uniform  distribution
sigma <- 1 # sigma > 0
x <- sqrt(-2*sigma^2*log(1-u)) # F(x) = 1 - exp(-x^2/(2*sigma^2)) , x >= 0
hist(x, prob = TRUE, main = expression(f(x) == (x / sigma^2)  * exp( - x^{2} / (2 * sigma ^ 2)) , sigma = 1),ylim = c(0,0.6))
y <- seq(0, 5, .01)
lines(y, (y / sigma^2) 
      * exp( - y^{2} / (2 * sigma ^ 2)))

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())
par(mfrow=c(1,2))
n <- 1000
u <- runif(n) # Uniform  distribution
sigma <- 2 # sigma > 0
x <- sqrt(-2*sigma^2*log(1-u)) # F(x) = 1 - exp(-x^2/(2*sigma^2)) , x >= 0
hist(x, prob = TRUE, main = expression(sigma == 2),ylim = c(0,0.4))
y <- seq(0, 8, .01)
lines(y, (y / sigma^2) 
      * exp( - y^{2} / (2 * sigma ^ 2)))

sigma <- 3 # sigma > 0
x <- sqrt(-2*sigma^2*log(1-u)) # F(x) = 1 - exp(-x^2/(2*sigma^2)) , x >= 0
hist(x, prob = TRUE, main = expression(sigma == 3),ylim = c(0,0.3))
y <- seq(0, 12, .01)
lines(y, (y / sigma^2) 
      * exp( - y^{2} / (2 * sigma ^ 2)))
layout(1)

## -----------------------------------------------------------------------------
rm(list=ls())
n <- 1000
X1 <- rnorm(n,0,sd=1)
X2 <- 3 + rnorm(n,0,sd=1)
p1 <- 0.75 # 0.50时有明显双峰
p2 <- 1 - p1
choice <- rbinom(n, 1, p1)
Z <- choice * X1 + (1-choice) * X2
par(mfrow=c(1,3))
hist(X1);hist(Z);hist(X2)
layout(1)

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())
n <- 1000
X1 <- rnorm(n,0,sd=1)
X2 <- 3 + rnorm(n,0,sd=1)
p1 <- 0.50 
p2 <- 1 - p1
choice <- rbinom(n, 1, p1)
Z <- choice * X1 + (1-choice) * X2
par(mfrow=c(1,3))
hist(X1);hist(Z);hist(X2)
layout(1)

## -----------------------------------------------------------------------------
rm(list=ls())
n <- 1000
t <- 10
lambda <- 1
r <- 1
beta <- 1
Nt <- rpois(n, lambda * t)
X <- rep(0,n)
for (i in 1:n) {
  Y <- rgamma(Nt[i], r, beta)
  X[i] <- sum(Y)
} # Now, X is the vector whose elements are the samples we need to use
tmean <- lambda * t * (r / beta)
tvariance <- lambda * t * ((r+1)*r / beta)

mean(X); tmean # tmean is the theoretical mean
var(X); tvariance # tvariance is the theoretical variance

## -----------------------------------------------------------------------------
rm(list=ls())
n <- 1000
X = 0.1*(1:10)
eY <- rep(0,10)
B = beta(3,3)
f <- function(x){
  (1/B) * x^2 * (1 - x)^2
}
for (i in 1:10) {
  u <- runif(n,0,X[i])
  g <- X[i] * f(u)
  eY[i] <- mean(g) #蒙特卡洛方法计算对应值
}
tY <- pbeta(X,3,3) #理论Beta分布函数对应值
par(mfrow = c(1,2))
plot(X,tY,type = "o",main = "理论值")
plot(X,eY,type = "o",main = "蒙特卡洛计算值")
layout(1)

## -----------------------------------------------------------------------------
rm(list=ls())
upper_value <- 5 #积分上限
sigma <- 1 #参数取值
Rayleigh <- function(x, R = 1000, antithetic = FALSE) {
  u <- runif(R/2)
  if (antithetic) v <- 1 - u else v <- runif(R/2)
  u <- c(u, v)
  g <- ((x*u) / sigma^2) * exp(-(x*u)^2 / (2*sigma^2)) # x*u ~ N(0,x)
  cdf <- mean(g)
  cdf
} # 1000个取值取均值
m <- 1000 # 1000次实验
Rayleigh1 <- Rayleigh2 <- numeric(m)
x <- 5
for (i in 1:m) {
  Rayleigh1[i] <- Rayleigh(x, R = 1000, antithetic = FALSE)#第一类，未对偶
  Rayleigh2[i] <- Rayleigh(x, R = 1000, antithetic = TRUE)#第二类，对偶
}
round(c(mean(Rayleigh1),mean(Rayleigh2)),4)
round(c(sd(Rayleigh1),sd(Rayleigh2),sd(Rayleigh2)/sd(Rayleigh1)),4)


## -----------------------------------------------------------------------------
# 调整方差的方法: importance function.
rm(list=ls())
m <- 1000 #调整方差后的蒙特卡洛方法实验次数
MC_times <- 1000 #每次实验的样本量
cdf0 <- cdf1 <- cdf2 <- numeric(MC_times)
theta.hat <- se <- numeric(3) # f0 = 1, f1, f2.

g <- function(x) {
  (1/(sqrt(2*pi)*x^4)) * exp(-1/(2*x^2)) * (x >= 0) * (x <= 1)
}

f2 <- function(x){
  (1/sqrt(2*pi)) *exp(- x^2 / 2)
} # 标准正态分布的p.d.f.

for (i in 1:MC_times) {
  x <- runif(m) #using f0
  fg <- g(x)
  cdf0[i] <- mean(fg) 
}
theta.hat[1] <- mean(cdf0) # f0下实验结果均值
se[1] <- sd(cdf0)# f0下实验结果方差

for (i in 1:MC_times) {
  x <- rexp(m, 1) #using f1
  fg <- g(x) / exp(-x)
  cdf1[i] <- mean(fg)
}
theta.hat[2] <- mean(cdf1)# f1下实验结果均值
se[2] <- sd(cdf1)# f1下实验结果方差

for (i in 1:MC_times) {
  x <- rnorm(n = m,mean = 0,sd = 1) #using f2
  fg <- g(x) / f2(x)
  cdf2[i] <- mean(fg)
}
theta.hat[3] <- mean(cdf2)# f1下实验结果均值
se[3] <- sd(cdf2)# f2下实验结果方差

round(theta.hat,4)

X = (1:1000)*.001
Y = g(X)
Y0 = rep(1,length(X))
Y1 = exp(-X)
Y2 = f2(X)
plot(X,Y,type = "l",ylim = c(0,1.4),main = "f(x)与g(x)的函数图像",col = 1)
lines(X,Y0,type = "l",col = 2)
lines(X,Y1,type = "l",col = 3)
lines(X,Y2,type = "l",col = 4)
legend(x = "topleft",
       legend = c("g(x)","f(x)=1","f(x)为参数为1指数函数",
                  "f(x)为标准正态分布"),
       col = c(1,2,3,4),lty = rep(1,4), cex = 0.6)

## -----------------------------------------------------------------------------
n = 20 # 总共可对n*1e4数目的数据排序
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}#??????
} # 快速排序法
An = rep(0,n)
Tn = An
for (i in 1:n) {
    test<-sample(x = 1:(n*1e4), size = i*1e4)
    An[i] <- system.time(quick_sort(test))[1] # 快速排序法所需时间
    Tn[i] <- (i*1e4)*log((i*1e4)) # nlogn
}
regression_result <- lm(An~Tn)
summary(regression_result)
plot(Tn,An,main = "scatter plot and regression line",xlab = "nlog(n)",ylab = "computation time")
En <- regression_result$coefficients[1] + regression_result$coefficients[2] * Tn
lines(Tn,En) # 回归线

## -----------------------------------------------------------------------------
rm(list=ls())
n_skewness = 1000; n_MC = 10000
# 生成数据的函数
generate_data <- function(n_skewness = 1000, n_MC = 10000){
  skewness <- numeric(n_skewness)
  for (i in 1:n_skewness) {
    e <- rnorm(n_MC)
    skewness[i] <- mean(((e-mean(e))/sd(e))^3)
  }
  return(skewness)
}
e <- generate_data(n_skewness,n_MC)
qt = c(0.025,0.05,0.95,0.975)
#求解分位数
Quantiles <- function(unknown_data, qt){
  quantile(unknown_data,qt)
}
query_x <- Quantiles(e,qt)
round(query_x,4)

## -----------------------------------------------------------------------------
SEE <- function(unknown_data,query_x,n_skewness){
  # 使用核密度估计估计密度函数
  density_est <- density(unknown_data)
  # 在密度估计的自变量值x中查找最接近的值
  L <- length(query_x)
  nearest_index <- numeric(L)
  for (i in 1:L) {
    nearest_index[i] <- which.min(abs(density_est$x - query_x[i]))
  }
  # 获取对应自变量值的密度值
  corresponding_density <- density_est$y[nearest_index]
  # 输出估计值的标准差
  SEE_result <- numeric(L)
  for (j in 1:L) {
    VAR_q <- qt[j]*(1-qt[j]) / (n_skewness * corresponding_density[j]^2)
    SEE_result[j] <- sqrt(VAR_q)
  }
  return(SEE_result)
}
SEE_result <- SEE(e,query_x,n_skewness)
round(SEE_result,4)

## ----eval=FALSE---------------------------------------------------------------
#  library(xtable)#引用xtable包
#  Comparation <- function(query_x , qt, n_MC){
#    # 渐进分布的分位数
#    a <- rnorm(n_MC,mean = 0,sd = sqrt(6/n_MC))
#    a_query_x <- quantile(a,qt)
#    Table <- rbind(query_x,a_query_x)
#    Table <- as.data.frame(Table)
#    colnames(Table) <- qt
#    row.names(Table) <- c("estimated quantiles",
#                          "the quantiles of the large sample approximation")
#    xtable(Table,digits = 4,caption = "Quantiles")
#  }
#  Comparation(query_x, qt, n_MC) #表格的latex代码

## -----------------------------------------------------------------------------
rm(list=ls())
library(MASS)

# set seed and create data vectors
set.seed(12345)
sample_size <- 1000 
sample_meanvector <- c(0, 0)
sample_covariance_matrix <- matrix(c(1,0.5,0.5,1),
                                   ncol = 2)
# create bivariate normal distribution
sample_distribution <- mvrnorm(n = sample_size,
                               mu = sample_meanvector, 
                               Sigma = sample_covariance_matrix)

## -----------------------------------------------------------------------------
alternative <- function(sample_size){
  u <- runif(sample_size)
  x <- u^{1/3} # F(x) = x^3, 0<=x<=1
  return(cbind(u,x))
}
X <- alternative(sample_size)

## -----------------------------------------------------------------------------
#correlation coefficient
Comparate_cor <- function(sample_distribution){
  COR <- numeric(3)
  COR[1] <- cor.test(sample_distribution[,1], sample_distribution[,2],
                method = "pearson")$estimate
  COR[2] <- cor.test(sample_distribution[,1], sample_distribution[,2],
                     method = "spearman")$estimate
  COR[3] <- cor.test(sample_distribution[,1], sample_distribution[,2],
                     method = "kendall")$estimate
  COR
}
# the result of experiments
# the sampled distribution is bivariate normal
COR1 <- Comparate_cor(sample_distribution)
# the sampled distribution is F(x) = x^3, 0<=x<=1
COR2 <- Comparate_cor(X)
# turn: "pearson" -> "spearman" -> "kendall"
round(COR1,4)
round(COR2,4)

## -----------------------------------------------------------------------------
rm(list=ls())
m <- 1e4
hypotheses_null_num <- 950
hypotheses_alternative_num <- 50
N <- hypotheses_null_num + hypotheses_alternative_num
labels_hypotheses <- numeric(N)
labels_hypotheses[1:hypotheses_null_num] <- 
  rep(0,hypotheses_null_num)
labels_hypotheses[(hypotheses_null_num+1):N] <-
  rep(1,hypotheses_alternative_num)
data_generate <- function(n1,n2,m){
  Data <- matrix(0, nrow = (n1+n2), ncol = m)
  for (i in 1:m) {
    data_null <- runif(n1)
    data_alternative <- rbeta(n2,0.1,1)
    data <- numeric(n1+n2)
    data[1:n1] <- data_null
    data[(n1+1):(n1+n2)] <- data_alternative
    Data[,i] <- data
  }
  return(Data)
}
Data <- data_generate(hypotheses_null_num,
                      hypotheses_alternative_num,m)
Bonferroni <- function(data,m,N){
  Data <- matrix(0, nrow = N, ncol = m)
  for (i in 1:m) {
    Data[,i] <- p.adjust(data[,i],method = 'bonferroni')
  }
  return(Data)
}
BH <- function(data,m,N){
  Data <- matrix(0, nrow = N, ncol = m)
  for (i in 1:m) {
    Data[,i] <- p.adjust(data[,i],method = 'BH')
  }
  return(Data)
}
alpha <- 0.1

Data.Bonferroni <- Bonferroni(Data,m,N)
Data.BH <- BH(Data,m,N)
FDR.Bonferroni <- FDR.BH <- TPR.Bonferroni <- TPR.BH <-
  FWER.Bonferroni <- FWER.BH <- numeric(m)
for (i in 1:m) {
  R <- sum(Data.Bonferroni[,i] <= alpha)
  V <- sum(Data.Bonferroni[1:hypotheses_null_num,i] <= alpha)
  S <- R-V
  U <- sum(Data.Bonferroni[(1+hypotheses_null_num):N,i] <= alpha)
  TT <- N-R-U
  FDR.Bonferroni[i] <- V/R
  TPR.Bonferroni[i] <- S/hypotheses_alternative_num
  if (V>0) {
    FWER.Bonferroni[i] <- 1
  } else {
    FWER.Bonferroni[i] <- 0
  }
}

for (i in 1:m) {
  R <- sum(Data.BH[,i] <= alpha)
  V <- sum(Data.BH[1:hypotheses_null_num,i] <= alpha)
  S <- R-V
  U <- sum(Data.BH[(1+hypotheses_null_num):N,i] <= alpha)
  TT <- N-R-U
  FDR.BH[i] <- V/R
  TPR.BH[i] <- S/hypotheses_alternative_num
  if (V>0) {
    FWER.BH[i] <- 1
  } else {
    FWER.BH[i] <- 0
  }
}

result <- matrix(0,nrow = 2,ncol = 3)
result[1,] <- c(mean(FWER.Bonferroni),mean(FDR.Bonferroni),mean(TPR.Bonferroni))
result[2,] <- c(mean(FWER.BH),mean(FDR.BH),mean(TPR.BH))
result <- as.data.frame(result,row.names = c("Bonferroni","BH"))
colnames(result) <- c("FWER","FDR","TPR")
result

## -----------------------------------------------------------------------------
library(boot)
rm(list=ls())
data <- as.array(aircondit$hours)
MLE <- function(x){
  return(1/mean(x))
} 
SIZE <- length(data)
R <- 1e4
Mboot <- replicate(R, expr = {
  y <- sample(data, size = SIZE, replace = TRUE)
  MLE(y)
})
round(c(bias = mean(Mboot) - MLE(data),se = sd(Mboot)),3)

## -----------------------------------------------------------------------------
library(boot)
rm(list=ls())
data <- as.array(aircondit$hours)
SIZE <- length(data)
m <- 1e2
boot.MLE <- function(x,i) mean(x[i])
ci.norm<-ci.basic<-ci.perc<-ci.bca<-matrix(NA,m,2)
for(i in 1:m){
  de <- boot(data,statistic=boot.MLE, R = 999)
  ci <- boot.ci(de,conf = 0.95,
                type=c("norm","basic","perc","bca")) # level: 0.95
  ci.norm[i,]<-ci$norm[2:3];ci.basic[i,]<-ci$basic[4:5]
  ci.perc[i,]<-ci$percent[4:5];ci.bca[i,]<-ci$bca[4:5]
}
# The standard bootstrap CI based on asymptotic normality
round(c(mean(ci.norm[,1]),mean(ci.norm[,2])),4) 
round(c(mean(ci.basic[,1]),mean(ci.basic[,2])),4) # The basic bootstrap CI
round(c(mean(ci.perc[,1]),mean(ci.perc[,2])),4) # Percentile CI
round(c(mean(ci.bca[,1]),mean(ci.bca[,2])),4) # Bias-corrected and accelerated CI

## -----------------------------------------------------------------------------
rm(list=ls())
library(bootstrap)
attach(scor)
head(scor)
B <- 1e3; set.seed(12345)
cor_matrix <- cor(scor)
straight_T <- function(cor_matrix){
  ev <- eigen(cor_matrix)
  max(ev$val[1]) / sum(ev$val)
} # 直接求出的参数估计

bootstrap_data <- function(data){
  x <- 1:nrow(data)
  xstar <- sample(x,replace=TRUE)
  ndata <- matrix(rep(0,nrow(data)*ncol(data)),nrow = nrow(data),ncol = ncol(data))
  for (i in 1:nrow(data)) ndata[i,] <- as.numeric(data[xstar[i],])
  ndata <- as.data.frame(ndata)
  colnames(ndata) <- colnames(data)
  return(ndata)
} # 重排数据生成
bootstrap_T <- function(data){
  total <- numeric(B)
  for(b in 1:B){
    ndata <- bootstrap_data(data)
    cor_matrix <- cor(ndata)
    ev <- eigen(cor_matrix)
    total[b] <- max(ev$val[1]) / sum(ev$val)
  }
  total
} # bootstrap方法求出的theta估计

jackknife_T <- function(data){
  total <- numeric(nrow(data))
  for (b in 1:nrow(data)) {
    ndata <- data[(1:nrow(data))[-b],]
    cor_matrix <- cor(ndata)
    ev <- eigen(cor_matrix)
    total[b] <- max(ev$val[1]) / sum(ev$val)
  }
  total
} # jackknife方法求出的theta估计

n <- nrow(scor)
theta <- straight_T(cor_matrix)
theta_bootstrap <- bootstrap_T(scor)
theta_jackknife <- jackknife_T(scor)
mean_jackknife <- rep(mean(theta_jackknife),nrow(scor))
se_jackknife <- (n-1)*mean((theta_jackknife-mean_jackknife)^2)
round(c(original=theta,bias.bootstrap=mean(theta_bootstrap)-theta,
        se.bootstrap=sd(theta_bootstrap),
        bias.jackknife=(n-1)*(mean(theta_jackknife)-theta),
        se.jackknife=se_jackknife),4)
detach(scor)

## -----------------------------------------------------------------------------
rm(list=ls())
library(DAAG); attach(ironslag)
a <- seq(10, 40, .1) #sequence for plotting fits

# 判据：调整的R方
L1 <- lm(magnetic ~ chemical) # 记为模型1
L2 <- lm(magnetic ~ chemical + I(chemical^2))# 记为模型2
L3 <- lm(magnetic ~ chemical + I(chemical^2) + I(chemical^3))# 记为模型3
L4 <- lm(log(magnetic) ~ log(chemical))# 记为模型4
summary1 <- summary(L1)
summary2 <- summary(L2)
summary3 <- summary(L3)
summary4 <- summary(L4)
ADJ.R.squared_result <- c(summary1$adj.r.squared,summary2$adj.r.squared,
               summary3$adj.r.squared,summary4$adj.r.squared)
which.min(ADJ.R.squared_result)

# for n-fold cross validation
# fit models on leave-one-out samples
CV <- function(M,C){
  n <- length(M) #in DAAG ironslag
  e1 <- e2 <- e3 <- e4 <- numeric(n)
  for (k in 1:n) {
    y <- M[-k]
    x <- C[-k]
    J1 <- lm(y ~ x)# 记为模型1
    yhat1 <- J1$coef[1] + J1$coef[2] * C[k]
    e1[k] <- M[k] - yhat1
    J2 <- lm(y ~ x + I(x^2))# 记为模型2
    yhat2 <- J2$coef[1] + J2$coef[2] * C[k] +
      J2$coef[3] * C[k]^2
    e2[k] <- M[k] - yhat2
    J3 <- lm(y ~ x + I(x^2) + I(x^3))# 记为模型3
    logyhat3 <- J3$coef[1] + J3$coef[2] * C[k]
    yhat3 <- exp(logyhat3)
    e3[k] <- M[k] - yhat3
    J4 <- lm(log(y) ~ log(x))# 记为模型4
    logyhat4 <- J4$coef[1] + J4$coef[2] * log(C[k])
    yhat4 <- exp(logyhat4)
    e4[k] <- M[k] - yhat4
  }
  result <- c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2))
  result
}

CV_result <- CV(M = magnetic, C = chemical)
which.min(CV_result)

## -----------------------------------------------------------------------------
rm(list=ls())
library(twosamples)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))

R <- 999 #number of replicates
z <- c(x, y) #pooled sample
K <- 1:26
reps <- numeric(R) #storage for replicates
t0 <- cvm_test(x,y)[1]
for (i in 1:R) {
  #generate indices k for the first sample
  k <- sample(K, size = 14, replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k] #complement of x1
  reps[i] <- cvm_test(x1, y1)[1]
}
p <- mean(c(t0, reps) <= t0)
detach(chickwts)
round(p,4)

## -----------------------------------------------------------------------------
rm(list=ls())
data <- as.data.frame(state.x77)
x <- data$Income
y <- data$Illiteracy
R <- 999 #number of replicates
n <- length(x)
reps <- numeric(R) #storage for replicates
result0 <- cor.test(x,y,method = "spearman",exact=FALSE)
t0 <- result0$p.value

for (i in 1:R) {
  #generate indices k for the first sample
  y1 <- sample(y, size = n, replace = FALSE) #complement of x1
  result <- cor.test(x,y1,method = "spearman",exact=FALSE)
  reps[i] <- result$p.value
}
p <- mean(c(t0, reps) >= t0)
round(p,4)

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(12345)
f <- function(x, theta,eta) {
  num <- theta * pi * (1 + ((x - eta)/theta)^2)
  1/num
} # theta > 0
m <- 10000
theta <- 1
eta <- 0
# 提议分布为N(X,sigma^2)
sigma <- 4
x <- numeric(m)
x[1] <- rnorm(n = 1,mean = 0,sd = sigma) 
k <- 0
u <- runif(m)
for (i in 2:m) {
  xt <- x[i-1]
  y <- rnorm(n = 1,mean = xt,sd = sigma)
  num <- f(y, theta, eta) * dnorm(xt,mean = y,sd = sigma)
  den <- f(xt, theta, eta) * dnorm(y,mean = xt,sd = sigma)
  if (u[i] <= num/den){
    x[i] <- y
  } else {
    x[i] <- xt
    k <- k+1     #y is rejected
  }
}
b <- 1001
y <- x[b:m]
Q <- quantile(y, probs=seq(0.1,1, by=0.1))
QC <- qt(df = 1,p = seq(0.1,1, by=0.1))
round(Q,3)
round(QC,3)

## -----------------------------------------------------------------------------
rm(list=ls())
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
a <- 1
b <- 1
n <- 10
###### generate the chain #####

X[1, ] <- c(n/2, 0.5) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
X[i, 1] <- rbinom(1,n,x2)
x1 <- X[i, 1]
X[i, 2] <- rbeta(1,x1+a,n-x1+b)
} # 已生成马尔科夫链

B <- burn + 1
x <- X[B:N, ]
cat('Means: ',round(colMeans(x),2))

plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c(expression(x),expression(y)),col=1:2,lwd=2)

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())
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
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}
f <- function(x, theta,eta) {
  num <- theta * pi * (1 + ((x - eta)/theta)^2)
  1/num
} # theta > 0
m <- 30000
theta <- 1
eta <- 0
# 提议分布为N(X,sigma^2)
sigma <- 4
norm_chain <- function(m,x0,sigma){
  x <- numeric(m)
  x[1] <- rnorm(n = 1,mean = x0,sd = sigma) 
  k <- 0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i-1]
    y <- rnorm(n = 1,mean = xt,sd = sigma)
    num <- f(y, theta, eta) * dnorm(xt,mean = y,sd = sigma)
    den <- f(xt, theta, eta) * dnorm(y,mean = xt,sd = sigma)
    if (u[i] <= num/den){
      x[i] <- y
    } else {
      x[i] <- xt
      k <- k+1     #y is rejected
    }
  }
  return(x)
}
k <- 4        #number of chains to generate
x0 <- c(-20,-10,10,20)
#generate the chains
#set.seed(12345)
X <- matrix(0, nrow=k, ncol=m)
for (i in 1:k)
  X[i, ] <- norm_chain(m,x0[i],sigma)

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

b <- 1000
for (i in 1:k){
  if(i==1){
    plot((b+1):m,psi[i, (b+1):m], type="l",
            xlab='Index', ylab=bquote(phi),ylim = c(-1,1))
  }else{
    lines(psi[i, (b+1):m], col=i)
  }
}
rhat <- rep(0, m)
for (j in (b+1):m)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):m], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## ----echo=FALSE---------------------------------------------------------------
rm(list=ls())
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
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}
m <- 30000
a <- 1
b <- 1
n <- 10
biviarate_chain <- function(N,x0){
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  X[1, ] <- x0 #initialize
  for (i in 2:N) {
    x2 <- X[i-1, 2]
    X[i, 1] <- rbinom(1,n,x2)
    x1 <- X[i, 1]
    X[i, 2] <- rbeta(1,x1+a,n-x1+b)
  } # 已生成马尔科夫链
  X
}
k <- 4        #number of chains to generate
x0 <- matrix(0, nrow=k, ncol=2)
x0[,1] <- runif(4,0,10)
#generate the chains
X <- matrix(0, nrow=k, ncol=m)
M <- matrix(0, nrow=m, ncol=2)
for (i in 1:k){
  M <- biviarate_chain(m,x0[i,])
  X[i, ] <- M[,1]
}

x0 <- matrix(0, nrow=k, ncol=2)
x0[,1] <- rep(n/2,4)
x0[,2] <- runif(4,0,1)
#generate the chains
Y <- matrix(0, nrow=k, ncol=m)
M <- matrix(0, nrow=m, ncol=2)
for (i in 1:k){
  M <- biviarate_chain(m,x0[i,])
  Y[i, ] <- M[,2]
}

#compute diagnostic statistics
psi1 <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi1))
  psi1[i,] <- psi1[i,] / (1:ncol(psi1))

b <- 1000
for (i in 1:k){
  if(i==1){
    plot((b+1):m,psi1[i, (b+1):m], type="l",
         xlab='Index', ylab=bquote(phi),ylim = c(4,6))
  }else{
    lines(psi1[i, (b+1):m], col=i)
  }
}
rhat <- rep(0, m)
for (j in (b+1):m)
  rhat[j] <- Gelman.Rubin(psi1[,1:j])
plot(rhat[(b+1):m], type="l", xlab="", ylab="R"
     ,main = "Y初始值不变，X初始值变化")
abline(h=1.2, lty=2)

#compute diagnostic statistics
psi2 <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psi2))
  psi2[i,] <- psi2[i,] / (1:ncol(psi2))

for (i in 1:k){
  if(i==1){
    plot((b+1):m,psi2[i, (b+1):m], type="l",
         xlab='Index', ylab=bquote(phi),ylim = c(0.4,0.6))
  }else{
    lines(psi2[i, (b+1):m], col=i)
  }
}
rhat <- rep(0, m)
for (j in (b+1):m)
  rhat[j] <- Gelman.Rubin(psi2[,1:j])
plot(rhat[(b+1):m], type="l", xlab="", ylab="R"
     ,main = "X初始值不变，Y初始值变化")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
rm(list=ls())
# 计算出项数为k_term，维数为dim，向量为vector_a的值
numerical <- function(k_term,dim,vector_a){
  factorial_k <- numeric(k_term)
  for (i in 1:k_term) factorial_k[i] <- log(i)
  # 先对数化
  term_1 <- -(k_term*log(2)+sum(factorial_k))
  norm_a <- sum(vector_a^2)
  term_2 <- (k_term+1)*log(norm_a)
  term_3 <- lgamma((dim+1)/2)+lgamma(k_term+1.5)-lgamma(k_term+dim/2+1)
  term <- term_1+term_2+term_3
  # 再指数化代入
  result <- (-1)^k_term * exp(term) / ((2*k_term+1)*(2*k_term+2))
  return(result)
}
dim <- 100; vector_a <- rep(2,dim); k_term <- 100
result <- numerical(k_term,dim,vector_a)
round(result,4)
# 项数求和，最终定会收敛
sum_numerical <- function(d,a){
  k <- 1
  minus <- numerical(k,d,a)
  Sum <- 0
  while (abs(minus) >= 1e-8) {
    Sum <- Sum + minus; k <- k+1 
    minus <- numerical(k,d,a)
  }
  Sum
}
d <- 2 ; a <- c(1,2); SUM <- sum_numerical(d,a)
round(SUM,4)


## -----------------------------------------------------------------------------
rm(list = ls())
f1 <- function(u, k) {
  below <- 1 + u ^ 2 / (k - 1)
  re <- -(k / 2) * exp(below)
  re
}
f2 <- function(u, k) {
  below <- 1 + u ^ 2 / k
  re <- -((k + 1) / 2) * exp(below)
  re
}
c <- function(a, k)
  sqrt(a ^ 2 * k / (k + 1 - a ^ 2))

rootfinding <- function(a, k) {
  coef <- sqrt(k / (k - 1)) * exp(2 * lgamma(k / 2) - lgamma((k - 1) / 2) -
                                    lgamma((k + 1) / 2))
  res1 <- integrate(f1,
                    lower = 0,
                    upper = c(a, k - 1),
                    k = k)
  res2 <- integrate(f2,
                    lower = 0,
                    upper = c(a, k),
                    k = k)
  re <- coef * res1$value - res2$value
  re
}
# 11.5结果：
K <- 25
res1 <- function(k) {
  i = 1
  # yleft <- 1; yright <- -1;
  xleft <- 0
  xright <- sqrt(k)
  mid <- rootfinding((xleft + xright) / 2, k)
  epsilon1 <- 1e-6
  de <- mid
  while (abs(de) >= epsilon1) {
    if (mid > 0) {
      xleft <- (xleft + xright) / 2
      de <- mid - rootfinding((xleft + xright) / 2, k)
      mid <- rootfinding((xleft + xright) / 2, k)
    }
    else {
      xright <- (xleft + xright) / 2
      de <- mid - rootfinding((xleft + xright) / 2, k)
      mid <- rootfinding((xleft + xright) / 2, k)
    }
    i <- i + 1
  }
  return((xleft + xright) / 2)
}
# 11.5的A(k):
root <- res1(K)
round(root, 4)

curves_intersection <- function(a, k) {
  S_lower <- 1 - pt(q = c(a, k - 1), df = k - 1)
  S_upper <- 1 - pt(q = c(a, k), df = k)
  return(S_lower - S_upper)
}
res2 <- function(k) {
  i = 1
  # yleft <- 1; yright <- -1;
  xleft <- 0
  xright <- sqrt(k)
  mid <- curves_intersection((xleft + xright) / 2, k)
  epsilon1 <- 1e-6
  de <- mid
  while (abs(de) >= epsilon1) {
    if (mid > 0) {
      xleft <- (xleft + xright) / 2
      de <- mid - curves_intersection((xleft + xright) / 2, k)
      mid <- curves_intersection((xleft + xright) / 2, k)
    }
    else {
      xright <- (xleft + xright) / 2
      de <- mid - curves_intersection((xleft + xright) / 2, k)
      mid <- curves_intersection((xleft + xright) / 2, k)
    }
    i <- i + 1
  }
  return((xleft + xright) / 2)
}
# 11.4的A(k):
root <- res2(K)
round(root, 4)


## -----------------------------------------------------------------------------
rm(list = ls())
y <- c(0.54,0.48,0.33,0.43,1.00,1.00,0.91,1.00,0.21,0.85)
tau <- 1
trim <- rep(1,length(y))
trim[which(y == 1)] <- 0
mlogL <- function(lambda = 1) {
  # minus log-likelihood
  n <- length(trim)
  t <- y
  t[which(trim == 0)] <- 0
  # 统计未被截断的观察值个数
  trimmed <- length(which(trim == 1))
  lambda <- abs(lambda)
  # cons + algebra 为似然函数
  cons <- - trimmed * log(lambda) - (n-trimmed) * (tau/lambda)
  algebra <- -sum(t)/lambda
  return(-(cons + algebra))
}
library(stats4)
fit <- mle(mlogL) # No additional data are allowed!
cat("MLE下的参数lambda估计值:",abs(as.numeric(fit@coef)))

# 初始化：
lambda0 <- 1

# EM算法实现
for (step in 1:200) {
  # E-step：
  op <- function(lambda){
    n <- length(trim)
    t <- y
    t[which(trim == 0)] <- 0
    trimmed <- length(which(trim == 1))
    lambda <- abs(lambda)
    # 计算条件期望
    delta <- choose(n-1,trimmed-1)*(1-exp(-tau/lambda0))^(trimmed-1)*(exp(-tau/lambda0))^(n-trimmed)
    delta <- delta * (1 - exp(-tau/lambda0))
    delta <- delta / (choose(n,trimmed)*(1-exp(-tau/lambda0))^(trimmed)*(exp(-tau/lambda0))^(n-trimmed))
    # fterm+sterm 是似然函数取条件期望的结果
    fterm <- delta * (-log(lambda)-mean(y)/lambda)
    sterm <- -(1-delta) * (tau/lambda)
    return(fterm+sterm)
  }
  oldlambda <- lambda0
  # M-step：最大化
  res <- optimize(op,lower=1e-5,upper=1e5,maximum=TRUE)
  lambda0 <- abs(res$objective)
  # 设阈值：当上一步迭代得到的参数与下一步迭代得到的参数变化很小，即认为收敛
  threshold <- 1e-5
  if (sum(abs(lambda0 - oldlambda)) < threshold) break
}
cat("EM算法下的参数lambda估计值:",abs(lambda0))

## -----------------------------------------------------------------------------
rm(list = ls())
library(Rglpk)
obj <- c(4,2,9); dir <- c("<=","<=",">=",">=",">=")
mat <- matrix(c(2,1,1,0,0,1,-1,0,1,0,1,3,0,0,1),nrow = 5)
rhs <- c(2,3,0,0,0)
max <- FALSE
Rglpk_solve_LP(obj, mat, dir, rhs,  max = max)

## -----------------------------------------------------------------------------
rm(list = ls())
attach(mtcars)
formulas <- list(mpg ~ disp,mpg ~ I(1 / disp),
                 mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
lapply(formulas, lm)

## -----------------------------------------------------------------------------
rm(list = ls())
attach(mtcars)
bootstraps <- lapply(1:10, function(i) {
    rows <- sample(1:nrow(mtcars), rep = TRUE)
    mtcars[rows, ]
})
lapply(1:10,function(i) lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp))
library(purrr)
fit_models <- map(bootstraps, ~lm(mpg ~ disp, data = .))
fit_models

## -----------------------------------------------------------------------------
rm(list = ls())
attach(mtcars)
rsq <- function(mod) summary(mod)$r.squared
# 3
formulas <- list(mpg ~ disp,mpg ~ I(1 / disp),
                 mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
ex3 <- lapply(formulas, lm)
lapply(ex3, rsq)
# 4
bootstraps <- lapply(1:10, function(i) {
    rows <- sample(1:nrow(mtcars), rep = TRUE)
    mtcars[rows, ]
})
ex4 <- lapply(1:10,function(i) lm(bootstraps[[i]]$mpg~bootstraps[[i]]$disp))
lapply(ex4, rsq)

## -----------------------------------------------------------------------------
rm(list = ls())
trials <- replicate(100,t.test(rpois(10, 10), rpois(7, 10)),simplify = FALSE)
sapply(1:100,function(i) trials[[i]]$p.value)
sapply(trials, "[[", "p.value")

## -----------------------------------------------------------------------------
# 创建一个自定义函数 parallel_lapply
parallel_lapply <- function(input_list, fun, output_type, ...) {
  if(output_type == "vector") {
    result <- vapply(input_list, fun, FUN.VALUE = numeric(1), ...)
  } else if(output_type == "matrix") {
    result <- matrix(unlist(Map(fun, input_list)), nrow = length(input_list))
  } else {
    stop("Invalid output type. Please specify 'vector' or 'matrix'.")
  }
  
  return(result)
}

# 示例输入列表
input_list <- list(1:3, 4:6, 7:9)

# 调用 parallel_lapply 函数来计算每个向量的平均值并存储在向量中
output_vector1 <- parallel_lapply(input_list,mean, "vector")
output_vector2 <- parallel_lapply(input_list,mean, "matrix")

# 输出结果
output_vector1
output_vector2

# 比较parallel_lapply与lapply
lapply(input_list,mean)

## -----------------------------------------------------------------------------
fast_chisq_test <- function(x, y) {
  obs <- table(x, y)
  n <- sum(obs)
  
  # 计算期望频数
  exp_freq <- outer(rowSums(obs), colSums(obs)) / n
  
  # 计算卡方统计量
  chisq_stat <- sum((obs - exp_freq)^2 / exp_freq)
  
  return(chisq_stat)
}

# 示例数据
x <- c(1, 2, 1, 2, 1)
y <- c(1, 1, 2, 2, 1)

# 调用 fast_chisq_test 函数计算卡方统计量
chisq_stat <- fast_chisq_test(x, y)

# 输出结果
chisq_stat

## -----------------------------------------------------------------------------
fast_table <- function(x, y) {
  max_x <- max(x)
  max_y <- max(y)
  
  counts <- matrix(0, nrow = max_x, ncol = max_y)
  
  for (i in 1:length(x)) {
    counts[x[i], y[i]] <- counts[x[i], y[i]] + 1
  }
  
  return(counts)
}

fast_chisq_test <- function(x, y) {
  obs <- fast_table(x, y)
  n <- sum(obs)
  
  # 计算期望频数
  exp_freq <- outer(rowSums(obs), colSums(obs)) / n
  
  # 计算卡方统计量
  chisq_stat <- sum((obs - exp_freq)^2 / exp_freq)
  
  return(chisq_stat)
}

# 示例数据
x <- c(1, 2, 1, 2, 1)
y <- c(1, 1, 2, 2, 1)

# 调用 fast_chisq_test 函数计算卡方统计量
chisq_stat <- fast_chisq_test(x, y)

# 输出结果
chisq_stat

## ----eval=FALSE---------------------------------------------------------------
#  rm(list=ls())
#  #initialize constants and parameters
#  N <- 5000 #length of chain
#  burn <- 1000 #burn-in length
#  X <- matrix(0, N, 2) #the chain, a bivariate sample
#  a <- 1
#  b <- 1
#  n <- 10
#  ###### generate the chain #####
#  
#  X[1, ] <- c(n/2, 0.5) #initialize
#  for (i in 2:N) {
#    x2 <- X[i-1, 2]
#    X[i, 1] <- rbinom(1,n,x2)
#    x1 <- X[i, 1]
#    X[i, 2] <- rbeta(1,x1+a,n-x1+b)
#  } # 已生成马尔科夫链
#  
#  B <- burn + 1
#  x <- X[B:N, ]
#  cat('Means: ',round(colMeans(x),2))
#  
#  plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
#  lines(x[,2],col=2,lwd=2)
#  legend('bottomright',c(expression(x),expression(y)),col=1:2,lwd=2)

## ----eval=FALSE---------------------------------------------------------------
#  gibbsSamplingR <- function(N, thin) {
#    mat <- matrix(nrow = N, ncol = 2)
#    n = 10; a = 1; b = 1; x = n/2; y = 0.5
#    for (i in 1:N) {
#      for (j in 1:thin) {
#        x = rbinom(1, n, y)
#        y = rbeta(1, x + a,n - x + b)
#      }
#      mat[i, ] <- c(x, y)
#    }
#    mat
#  }

## -----------------------------------------------------------------------------
rm(list=ls())
library(Rcpp)
library(SA24204187)
N = 500; thin <- 10
GibbsR <- gibbsSamplingR(N,thin)
GibbsC <- gibbsSamplingC(N,thin)
par(mfrow=c(1,2))
plot(GibbsR[,1],type='l',col=1,lwd=2,xlab='Index-GibbsR',ylab='Random numbers')
lines(GibbsR[,2],col=2,lwd=2)
plot(GibbsC[,1],type='l',col=1,lwd=2,xlab='Index-GibbsR',ylab='Random numbers')
lines(GibbsC[,2],col=2,lwd=2)
layout(1)

## -----------------------------------------------------------------------------
rm(list=ls())
library(Rcpp)
library(SA24204187)
N = 500; thin <- 10
GibbsR <- gibbsSamplingR(N,thin)
GibbsC <- gibbsSamplingC(N,thin)
library(microbenchmark)
ts <- microbenchmark(GibbsR = gibbsSamplingR(N,thin),GibbsC = gibbsSamplingC(N,thin))
summary(ts)[,c(1,3,5,6)]

