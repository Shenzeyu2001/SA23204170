## -----------------------------------------------------------------------------
library(SA23204170)#之前作业用了一个不在R包里的数据集，所以这里要library一下自己的包
data("height_weight")
lm <- lm(weight~height,data = log(data))
summary(lm)$coef

## ----echo=FALSE---------------------------------------------------------------
lm <- lm(weight~height,data = log(data))
summary(lm)$coef

## -----------------------------------------------------------------------------
plot(lm)

## -----------------------------------------------------------------------------
my.sample<-function(data,size,prob=NULL){
  if(is.null(prob)){
    length <- length(data)
    prob <- (numeric(length)+1)/length
  }
  
  stopifnot(sum(prob)==1)#当给定的概率分布的和不为1时报错
  # 生成[0, 1)之间的均匀随机数
  U <- runif(size)
  # 计算累积概率
  cumulative_prob <- cumsum(prob)
  # 找到所在区间
  index <- findInterval(U,cumulative_prob)+1
  # 根据索引获取抽样的数据
  sampled_data <- data[index]
  
  return(sampled_data)
}

## -----------------------------------------------------------------------------
# 不给定概率分布
my.sample(letters,10)

# 给定概率分布
my.sample(letters,10,prob = c(rep(.1,5),rep(.5/21,21)))


## -----------------------------------------------------------------------------
n <- 1000
u <- runif(n)
index <- which(u >= 0.5)
x <- c(-log(2-2 * u[index]), log(2 * u[-index]))
hist(x, breaks = "Scott", prob = TRUE, ylim = c(0, 0.5),xlim = c(-7,7))
# 经验密度函数
lines(density(x),lwd=2,col="black")
# 理论密度函数
x_positive<-c(0, qexp(ppoints(100), rate = 1))
x_negative<--rev(x_positive)
lines(x_negative, 0.5 * exp(x_negative),lwd=2,col="red")
lines(x_positive, 0.5 * exp(-x_positive),lwd=2,col="red")

## -----------------------------------------------------------------------------
my.rbeta <- function(n, a, b) {
n <- 1000
k <- 0
y <- numeric(n)
while (k < n) {
  u <- runif(1)
  x <- runif(1)
    if (x^(a - 1) * (1 - x)^(b - 1) > u) {
      k <- k + 1
      y[k] <- x
    }
 }
return(y)
} 

## -----------------------------------------------------------------------------
beta_sample<-my.rbeta(1000,3,2)
hist(beta_sample, breaks = "Scott", prob = TRUE)
# 经验密度函数
lines(density(beta_sample),lwd=2,col="black")
# 理论密度函数
y <- seq(0, 1, 0.01)
fy <- 12 * y^2 * (1 - y)
lines(y, fy,lwd=2,col="red")

## -----------------------------------------------------------------------------
my.gen<-function(n){
  x<-rep(0,n)
  u1<-runif(n,-1,1)
  u2<-runif(n,-1,1)
  u3<-runif(n,-1,1)
  for (i in 1:n) {
    if ((abs(u3[i]) >= abs(u2[i])) && (abs(u3[i]) >= abs(u1[i])))
      x[i] <- u2[i]
    else x[i] <- u3[i]
  }
  return(x)
}

x<-my.gen(10000)
hist(x, breaks = "Scott", prob = TRUE)
# 经验密度函数
lines(density(x),lwd=2,col="black")
# 理论密度函数
y <- seq(-1, 1, 0.01)
fy <- 0.75 * (1 - y^2)
lines(y, fy,lwd=2,col="red")


## -----------------------------------------------------------------------------
# 把整个蒙特卡罗计算pi估计值的的过程整合为一个函数
var_pihat <- function(rho){
  n <- 10^6
  d <- 10
  l <- 10*rho
  pihat <- numeric(100)
  for (i in 1:100) {
    X <- runif(n,0,d/2)
    Y <- runif(n,0,pi/2)
    pihat[i] <- 2*l/d/mean(l/2*sin(Y)>X)
  }
  var(pihat)
}

## -----------------------------------------------------------------------------
# rho=0.3
var_pihat(0.3)
# rho=0.6
var_pihat(0.6)
# rho=1
var_pihat(1)

## -----------------------------------------------------------------------------
v_mc <- .5*(.5*(exp(2)-1)-(exp(1)-1)^2)
v_anti <-.5*(.5*(exp(2)-1)-(exp(1)-1)^2+exp(1)-(exp(1)-1)^2)
(v_mc-v_anti)/v_mc                                         

## -----------------------------------------------------------------------------
m <- 10000
mc <- replicate(1000, expr = {
  mean(exp(runif(m)))
  })
anti <- replicate(1000, expr = {
  U <- runif(m/2)
  mean((exp(U)+exp(1-U))/2)
  })
# 两种方法得到的估计值
c(mean(mc), mean(anti))
# 两种方法的方差
c(var(mc), var(anti))
# 方差减少的比例
(var(mc) - var(anti))/var(mc)

## -----------------------------------------------------------------------------
x <- seq(1, 10, .01)
g <- x^2 * exp(-x^2/2)/sqrt(2 * pi)
plot(x, g, ylab = "", type = "l", ylim = c(0, 1), lwd=2, col=1)
lines(x, 2*dnorm(x, 1, 1), lty = 2, lwd=2, col=2)
lines(x, 2*dt(x-1, 1), lty = 3, lwd=2,col=3)
legend("topright", inset = 0.02, legend = c("g(x)", "f1", "f2"), lty = 1:3, lwd=2, col=1:3)

## -----------------------------------------------------------------------------
plot(x,g/(2*dnorm(x, 1, 1)), ylab = "", type = "l", ylim = c(0, 1), lty = 2, lwd=2, col=2)
lines(x, g/(2*dt(x-1, 1)), lty = 3, lwd=2,col=3)
legend("topright", inset = 0.02, legend = c("f1", "f2"), lty = 2:3, lwd=2, col=2:3)

## -----------------------------------------------------------------------------
est <- sd <- numeric(2)
u <- runif(1e6)
g <- function(x){x^2 * exp(-x^2/2)/sqrt(2 * pi)}
x <- abs(rnorm(1e6))+1 #使用f1
fg <- g(x)/(2*dnorm(x, 1, 1))
est[1] <- mean(fg)
sd[1] <- sd(fg)
x <- abs(rt(1e6,1))+1
fg <- g(x)/(2*dt(x-1, 1))
est[2] <- mean(fg)
sd[2] <- sd(fg)
names(sd) <- c("f1","f2")
sd

## -----------------------------------------------------------------------------
est[1]

## -----------------------------------------------------------------------------
k <- 5
g <- function(x) exp(-x)/(1 + x^2)
f <- function(x){k/(1 - exp(-1))*exp(-x)}
M <- 1e5
m <- M/k
theta <- v <- numeric(5)
for (j in 1:k) {
  u <- runif(m, (j - 1)/k, j/k)
  x <- -log(1 - (1 - exp(-1)) * u)
  fg <- g(x)/f(x)
  theta[j] <- mean(fg)
  v[j] <- var(fg)
}
sum(theta)
sqrt(mean(v))

## -----------------------------------------------------------------------------
f <- function(x){1/(1 - exp(-1))*exp(-x)}
u <- runif(M)
x <- -log(1 - (1 - exp(-1)) * u)
fg <- g(x)/f(x)
theta <- mean(fg)
sd <- sd(fg)
theta
sd

## -----------------------------------------------------------------------------
n <- 20
t0 <- qt(c(0.025, 0.975), df = n - 1)
CI <- replicate(10000, expr = {
  x <- rchisq(n, df = 2)
  ci <- mean(x) + t0 * sd(x)/sqrt(n)
})
LCL <- CI[1,]
UCL <- CI[2,]
mean(LCL < 2 & UCL > 2)

## -----------------------------------------------------------------------------
test <- function(distribution){
  m <- 10000
  alpha <- 0.05
  count <- 0
  if(distribution=="chisq"){
    for (i in 1:m) {
     data <- rchisq(30, df = 1)
      result <- t.test(data, mu = 1, alternative = "two.sided")
      if (result$p.value < alpha) {
        count <- count + 1
      } 
    }
  }
  if(distribution=="unif"){
    for (i in 1:m) {
     data <- runif(30, min = 0, max = 2)
      result <- t.test(data, mu = 1, alternative = "two.sided")
      if (result$p.value < alpha) {
        count <- count + 1
      } 
    }
  }
  if(distribution=="exp"){
    for (i in 1:m) {
     data <- rexp(30, rate = 1) 
      result <- t.test(data, mu = 1, alternative = "two.sided")
      if (result$p.value < alpha) {
        count <- count + 1
      } 
    }
  }
  return(count/m)
}

## -----------------------------------------------------------------------------
test("chisq")

## -----------------------------------------------------------------------------
test("unif")

## -----------------------------------------------------------------------------
test("exp")

## -----------------------------------------------------------------------------
m <- 1000
alpha <- .1

rep<- replicate(1000,{
  H0 <- runif(.95*m,0,1)
  H1 <- rbeta(.05*m,.1,1)
  p <- c(H0,H1)
  p.adj1 <- p.adjust(p,method='BH')
  p.adj2 <- p.adjust(p,method='bonferroni')

  FWER <- FDR <- TPR <- numeric(2)
  FDR[1] <- sum(p.adj1[1:950]<alpha)/sum(p.adj1<alpha)
  FDR[2] <- sum(p.adj2[1:950]<alpha)/sum(p.adj2<alpha)
  FWER[1] <- (sum(p.adj1[1:950]<alpha)>0)
  FWER[2] <- (sum(p.adj2[1:950]<alpha)>0)
  TPR[1] <- sum(p.adj1[-(1:950)]<alpha)/50
  TPR[2] <- sum(p.adj2[-(1:950)]<alpha)/50
  c(FWER,FDR,TPR)
})
result <- apply(rep,1,FUN=mean)

A <- matrix(round(result,3),ncol=3,nrow = 2)
colnames(A)<- c("FWER","FDR","TPR")
rownames(A)<- c("B-H","Bonferroni")
knitr::kable(A, format = "latex")

## -----------------------------------------------------------------------------

lambda <- 2          # 真实参数 lambda
size <- c(5, 10, 20)  # 样本大小
B <- 1000                 # 重抽样次数
m <- 1000                 # boot次数

mean_bias <- numeric(length(size))      
mean_sd <- numeric(length(size))  

# 模拟循环
for (i in 1:length(size)) {
  n <- size[i]   # 当前样本大小
  bias <- numeric(m)   # 每次模拟的bias
  sd <- numeric(m)  # 每次模拟的sd
  
  for (j in 1:m) {
    # 生成样本
    data <- rexp(n, rate = 1/lambda)
    
    # Boot
    boot <- replicate(B, {
      data_b <- sample(data, replace = TRUE)
      lambda_b <- 1 / mean(data_b)
      lambda_b
    })
    
    # 计算偏差和标准误差
    lambda_hat <- 1/mean(data)
    bias[j] <- mean(boot) - lambda_hat
    sd[j] <- sd(boot)
  }
  
  mean_bias[i] <- mean(bias)
  mean_sd[i] <- mean(sd)
}

# 理论值
theo_bias <- lambda / (size - 1)
theo_sd <- lambda * sqrt(size) / ((size - 1) * sqrt(size - 2))


result <- matrix(round(c(size,mean_bias,theo_bias,mean_sd,theo_sd),3),nrow=3)
colnames(result) <- c("SIze","sim_bias","theo_bias","sim_sd","theo_sd")
knitr::kable(result, format = "latex")


## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
attach(law)  

cor.stat <- function(x,i) {
  cor(x[i, 1], x[i, 2])  # 计算相关性统计量
}

cor.stat2 <- function(x,i) {
  o <- boot(x[i, ], cor.stat, R = 100)
  n <- length(i)  # 样本大小
  # 返回相关性统计量和方差估计
  c(o$t0, var(o$t) * (n - 1) / n^2)
}

# 使用boot函数执行Bootstrap分析
b <- boot(law, statistic = cor.stat2, R = 1000)

result <- boot.ci(b, type = "stud")
result

## -----------------------------------------------------------------------------
library(boot)
x <- aircondit
meant <- function(x, i) return(mean(as.matrix(x[i, ])))
b <- boot(x, statistic = meant, R = 2000)
boot.ci(b, type = c("norm", "basic", "perc", "bca"))

## -----------------------------------------------------------------------------
library(bootstrap)
x <- as.matrix(scor)
n <- nrow(x)
theta.jack <- numeric(n)
lambda <- eigen(cov(x))$values
theta.hat <- max(lambda/sum(lambda))
for (i in 1:n) {
  y <- x[-i, ]
  s <- cov(y)
  lambda <- eigen(s)$values
  theta.jack[i] <- max(lambda/sum(lambda))
  }
bias.jack <- (n - 1) * (mean(theta.jack) - theta.hat)
se.jack <- sqrt((n - 1)/n * sum((theta.jack - mean(theta.jack))^2))
result <- c(theta.hat, bias.jack, se.jack)
A <- matrix(round(result,3),ncol=3,nrow = 1)
colnames(A)<- c("Estimate","Bias","Standard Error")
knitr::kable(A, format = "latex")

## -----------------------------------------------------------------------------
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)

R <- 999
n <- length(x)
m <- length(y)
z <- c(x,y)
N <- n + m
K <- Fn <- Gm <- numeric(N)
for (i in 1:N) {
  Fn[i] <- mean(as.integer(z[i] <= x))
  Gm[i] <- mean(as.integer(z[i] <= y))
}
stat0 <- ((n * m)/N) * sum((Fn - Gm)^2)
stat <- replicate(R,{
  K <- sample(1:N)
  Z <- z[K]
  X <- Z[1:n]
  Y <- Z[(n + 1):N]
  for (i in 1:N) {
    Fn[i] <- mean(as.integer(Z[i] <= X))
    Gm[i] <- mean(as.integer(Z[i] <= Y))
    }
  ((n * m)/N) * sum((Fn - Gm)^2)
})
stat1 <-c(stat0,stat)
cat("The statistic is",round(stat0,3),".","\n")
cat("The p-value is",round(mean(stat1 >= stat0),3),".")

## -----------------------------------------------------------------------------
maxout <- function(x,y){
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx,outy)))
}

max.text <- function(x,y,R=999){
  n <- length(x)
  z <- c(x,y)
  N <- length(z)
  stat0 <- maxout(x,y)
  stat <- replicate(R,{
    k <- sample(1:N)
    k1 <- k[1:n]
    k2 <- k[(n + 1):N]
    maxout(z[k1], z[k2])
  })
  stat1 <- c(stat0,stat)
  p_value <- mean(stat1 >= stat0)
  return(p_value)
}


## -----------------------------------------------------------------------------
# 方差相同
p <- replicate(1000,{
  x <- rnorm(20,1,1)
  y <- rnorm(40,1,1)
  max.text(x,y)
})
power<- sum(p<.05)/1000
cat("Power is",power,".\n")

## -----------------------------------------------------------------------------
# 方差不同
p <- replicate(1000,{
  x <- rnorm(20,1,1)
  y <- rnorm(40,1,2)
  max.text(x,y)
})
power<- sum(p<.05)/1000
cat("Power is",power,".\n")

## -----------------------------------------------------------------------------
alpha <- function(N,b1,b2,b3,f0){
  x1 <- rpois(N,1)
  x2 <- rexp(N)
  x3 <- rbinom(N,1,0.5)
  g <- function(alpha){
    tmp <- exp(-alpha-b1*x1-b2*x2-b3*x3)
    p <- 1/(1+tmp)
    mean(p) - f0
  }
  solution <- uniroot(g,c(-20,10))
  return(alpha=solution$root)
}

## -----------------------------------------------------------------------------
N <- 10^6
b1 <-0
b2 <- 1
b3 <- -1
f0 <- c(.1,.01,.001,.0001)
alpha(N,b1,b2,b3,f0[1])
alpha(N,b1,b2,b3,f0[2])
alpha(N,b1,b2,b3,f0[3])
alpha(N,b1,b2,b3,f0[4])

## -----------------------------------------------------------------------------
neg_logf0 <- c(seq(0.01,2,length.out=10),seq(2,10,length.out=20))
alpha_seq <- numeric(30)
for (i in 1:30) {
  alpha_seq[i] <- alpha(N,b1,b2,b3,exp(-neg_logf0[i]))
}
plot(neg_logf0,alpha_seq,xlab = expression(-logf_0),ylab = "alpha",type = "l")

## -----------------------------------------------------------------------------
set.seed(12345)
    rw.Metropolis <- function(sigma, x0, N) {
        x <- numeric(N)
        x[1] <- x0
        u <- runif(N)
        k <- 0
        for (i in 2:N) {
            y <- rnorm(1, x[i-1], sigma)
                if (u[i] <= exp(abs(x[i-1]) - abs(y))){
                   x[i] <- y  
                   k <- k + 1
                  }
                else {
                    x[i] <- x[i-1]
                }
            }
        return(list(x=x, k=k))
    }

    
    N <- 2000
    sigma <- c(.05, .5, 4,  16)

    x0 <- 25
    rw1 <- rw.Metropolis(sigma[1], x0, N)
    rw2 <- rw.Metropolis(sigma[2], x0, N)
    rw3 <- rw.Metropolis(sigma[3], x0, N)
    rw4 <- rw.Metropolis(sigma[4], x0, N)

    
    #number of candidate points accepted
    accept.rate <- data.frame(sigma=sigma,accept.rate=c(rw1$k, rw2$k, rw3$k, rw4$k)/N)
    knitr::kable(accept.rate,format='latex')

## -----------------------------------------------------------------------------
      #display 4 graphs together
    refline <- c(-qexp(.975),qexp(.975))
    rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
    for (j in 1:4) {
        plot(rw[,j], type="l",
             xlab=bquote(sigma == .(round(sigma[j],3))),
             ylab="X", ylim=range(rw[,j]))
        abline(h=refline)
    }
     #reset to default

## -----------------------------------------------------------------------------
  p <- ppoints(200)
  y <- qexp(p, 1)
  z <- c(-rev(y), y)
  b <- 200
  y1 <- rw1$x[(b + 1):N]
  y2 <- rw2$x[(b + 1):N]
  y3 <- rw3$x[(b + 1):N]
  y4 <- rw4$x[(b + 1):N]
  
  Q1 <- quantile(y1, p)
  qqplot(z, Q1, main=bquote(sigma == .(round(sigma[1],3))),
        xlab="Laplace Quantiles", ylab="Sample Quantiles")
  abline(0, 1,col='blue',lwd=2)
  Q2 <- quantile(y2, p)
  qqplot(z, Q2, main=bquote(sigma == .(round(sigma[2],3))),
        xlab="Laplace Quantiles", ylab="Sample Quantiles")
  abline(0, 1,col='blue',lwd=2)
  Q3 <- quantile(y3, p)
  qqplot(z, Q3, main=bquote(sigma == .(round(sigma[3],3))),
        xlab="Laplace Quantiles", ylab="Sample Quantiles")
  abline(0, 1,col='blue',lwd=2)
  Q4 <- quantile(y4, p)
  qqplot(z, Q4, main=bquote(sigma == .(round(sigma[4],3))),
        xlab="Laplace Quantiles", ylab="Sample Quantiles")
  abline(0, 1,col='blue',lwd=2)
  

## -----------------------------------------------------------------------------
#initialize constants and parameters
N <- 5000 #length of chain
burn <- 1000 #burn-in length
X <- matrix(0, N, 2) #the chain, a bivariate sample
rho <- .9 #correlation
mu1 <- 0
mu2 <- 0
sigma1 <- 1
sigma2 <- 1
s1 <- sqrt(1-rho^2)*sigma1
s2 <- sqrt(1-rho^2)*sigma2
###### generate the chain #####

X[1, ] <- c(mu1, mu2) #initialize
for (i in 2:N) {
x2 <- X[i-1, 2]
m1 <- mu1 + rho * (x2 - mu2) * sigma1/sigma2
X[i, 1] <- rnorm(1, m1, s1)
x1 <- X[i, 1]
m2 <- mu2 + rho * (x1 - mu1) * sigma2/sigma1
X[i, 2] <- rnorm(1, m2, s2)
}
b <- burn + 1
x <- X[b:N, ]
plot(x[,1],type='l',col=1,lwd=2,xlab='Index',ylab='Random numbers')
lines(x[,2],col=2,lwd=2)
legend('bottomright',c(expression(X),expression(Y)),col=1:2,lwd=2)

## -----------------------------------------------------------------------------
X <- x[, 1]
Y <- x[, 2]
fit <- lm(Y ~ X)
summary(fit)

## -----------------------------------------------------------------------------
qqnorm(fit$res)
qqline(fit$res, col=2,lwd=2)

## -----------------------------------------------------------------------------
plot(fit$fit, fit$res, cex = 0.5)
 abline(h = 0)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
# psi[i,j] is the statistic psi(X[i,1:j])
# for chain in i-th row of X
psi <- as.matrix(psi)
n <- ncol(psi)
k <- nrow(psi)
psi.means <- rowMeans(psi) #row means
B <- n * var(psi.means) #between variance est.
psi.w <- apply(psi, 1, "var") #within variances
W <- mean(psi.w) #within est.
v.hat <- W*(n-1)/n + (B/n) #upper variance est.
r.hat <- v.hat / W #G-R statistic
return(r.hat)
}

## -----------------------------------------------------------------------------
f <- function(x, sigma) {
if (x < 0)
return(0)
stopifnot(sigma > 0)
return((x/sigma^2) * exp(-x^2/(2 * sigma^2)))
}
Rayleigh.MH.chain <- function(sigma, m, x0) {
  x <- numeric(m)
  x[1] <- x0
  u <- runif(m)
  for (i in 2:m) {
    xt <- x[i - 1]
    y <- rchisq(1, df = xt)
    num <- f(y, sigma) * dchisq(xt, df = y)
    den <- f(xt, sigma) * dchisq(y, df = xt)
    if (u[i] <= num/den)
      x[i] <- y
    else x[i] <- xt
  }
  return(x)
}
sigma <- 4
x0 <- c(1/sigma^2, 1/sigma, sigma^2, sigma^3)
k <- 4
m <- 2000
X <- matrix(0, nrow = k, ncol = m)
for (i in 1:k) X[i, ] <- Rayleigh.MH.chain(sigma, m, x0[i])
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi)) psi[i, ] <- psi[i, ]/(1:ncol(psi))
rhat <- Gelman.Rubin(psi)
rhat

## -----------------------------------------------------------------------------
library(coda)
X1 <- as.mcmc(X[1, ])
X2 <- as.mcmc(X[2, ])
X3 <- as.mcmc(X[3, ])
X4 <- as.mcmc(X[4, ])
Y <- mcmc.list(X1, X2, X3, X4)
print(gelman.diag(Y))
gelman.plot(Y, col = c(1, 1))

## -----------------------------------------------------------------------------
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- u+1
# MLE 算法
L <- function(lambda){
  sum((v*exp(-lambda*v)-u*exp(-lambda*u))/(exp(-lambda*u)-exp(-lambda*v)))
}
lambda_mle <- uniroot(L,c(0,10))$root
cat("MLE算法得到的估计为",round(lambda_mle,4))

## -----------------------------------------------------------------------------
# E-M 算法
n <- length(u)
f <- function(lambda){
  sum(lambda*(u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v)))
}
g <- function(lambda){
  n*lambda/(f(lambda)+n)
}
lambda0 <- .000001 #初始值
iter <- 1
repeat{
  lambda_EM <- g(lambda0)
  iter <- iter+1
  if((abs(lambda_EM-lambda0)<.Machine$double.eps^0.25) || (iter > 1000)) break
  lambda0 <- lambda_EM
}
cat("E-M算法得到的估计为",round(lambda_EM,4))

## -----------------------------------------------------------------------------
solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}
A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
               2,0,0,0,-3,-3,4,0,0,
               2,0,0,3,0,0,0,-4,-4,
               -3,0,-3,0,4,0,0,5,0,
               0,3,0,-4,0,-4,0,5,0,
               0,3,0,0,4,0,-5,0,-5,
               -4,-4,0,0,0,5,0,0,6,
               0,0,4,-5,-5,0,0,0,6,
               0,0,4,0,0,5,-6,-6,0), 9, 9)
library(boot) #needed for simplex function
B <- A+2
s0 <- solve.game(A)
s <- solve.game(B)
s0$v # A的value
s$v # B的value
round(cbind(s$x, s$y),7)
round(s$x*61, 7)

## -----------------------------------------------------------------------------
X <-cars
dim(X[FALSE,])
dim(X[,FALSE])
dim(X[FALSE,FALSE])

## ----eval=FALSE---------------------------------------------------------------
#  scale01 <- function(x) {
#  rng <- range(x, na.rm = TRUE)
#  (x - rng[1]) / (rng[2] - rng[1])
#  }

## ----eval=FALSE---------------------------------------------------------------
#  data.frame(lapply(X, function(x) if (is.numeric(x)) scale01(x) else x))

## ----eval=FALSE---------------------------------------------------------------
#  vapply(X, sd, numeric(1))

## ----eval=FALSE---------------------------------------------------------------
#  vapply(X[vapply(X, is.numeric, logical(1))],
#         sd,
#         numeric(1))

## -----------------------------------------------------------------------------
N <- 10000
n <- 10
a <- 2
b <- 3
gibbsR <- function(N, n, a, b) {
  mat <- matrix(nrow = N, ncol = 2)
  x <- y <- 0
  for (i in 1:N) {
    x <- rbinom(1, n, y)
    y <- rbeta(1, x+a, n-x+b)
    mat[i, ] <- c(x, y)
  }
  mat
}

## -----------------------------------------------------------------------------
library(Rcpp)
dir_cpp <- '../HW9/'
# Can create source file in Rstudio
#sourceCpp(paste0(dir_cpp,"gibbsC.cpp"))
library(microbenchmark)
ts <- microbenchmark(gibbsR=gibbsR(N, n, a, b),gibbsC=gibbsC(N, n, a, b))
summary(ts)[,c(1,3,5,6)]

