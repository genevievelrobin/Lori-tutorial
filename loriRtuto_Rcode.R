library(lori)
set.seed(123)

## covariates
m1 <- 30 # number of rows
m2 <- 10 # number of columns
K1 <- 2 # number of row covariates
K2 <- 2 # number of column covariates
K3 <- 3 # number of (rowxcolumn) covariates
R <- matrix(rnorm(m1*K1), nrow=m1) # matrix of row covariates
C <- matrix(rnorm(m2*K2), nrow=m2) # matrix of column covariates
E <- matrix(rnorm(m1*m2*K3), nrow=m1*m2) # matrix of  (rowxcolumn) covariates
U <- covmat(m1, m2, R, C, E)
U <- scale(U)

## parameters
alpha0 <- rep(0, m1)
alpha0[1:6] <- 1
beta0 <- rep(0, m2)
beta0[1:4] <- 1
epsilon0 <- rep(0, K1+K2+K3)
epsilon0[5:6] <- 0.2
r <- 2 #rank of interaction matrix theta0
theta0 <- 0.1*matrix(rnorm(m1*r), nrow=m1)%*%diag(r:1)%*%matrix(rnorm(m2*r), nrow=r)
theta0 <- sweep(theta0, 2, colMeans(theta0))
theta0 <- sweep(theta0, 1, rowMeans(theta0))

## construct x0
x0 <- matrix(rep(alpha0, m2), nrow=m1) #row effects
x0 <- x0 + matrix(rep(beta0, each=m1), nrow=m1) #add column effects
x0 <- x0 + matrix(U%*% epsilon0, nrow=m1) #add cov effects
x0 <- x0 + theta0 #add interactions

## sample count data y
y0 <- matrix(rpois(m1*m2, lambda = c(exp(x0))), nrow = m1)

## add missing values
p <- 0.2
y <- y0
y[sample(1:(m1*m2), round(p*m1*m2))] <- NA

## lori estimation
lambda1 <- 0.1
lambda2 <- 0.1
m <- sum(!is.na(y))
t <- Sys.time()
res.lori <- lori(y, U, 0.1, 0.1, trace.it=T)
t <- Sys.time()-t


## cross-validation
res.cv <- cv.lori(y, U, trace.it = T)
res.lori <- lori(y, U, res.cv$lambda1, res.cv$lambda2)
#res.lori <- lori(y, U, res.cv$lambda1, res.cv$lambda2, reff=F)
#res.lori <- lori(y, U, res.cv$lambda1, res.cv$lambda2, ceff=F)
#res.lori <- lori(y, U, res.cv$lambda1, res.cv$lambda2, reff=F, ceff=F)

res.mi <- mi.lori(y, U, res.cv$lambda1, res.cv$lambda2)
res.pool <- pool.lori(res.mi)

boxplot(res.mi$mi.alpha, pch="")
par(cex.axis=1)
par(cex.lab=1)
boxplot(res.mi$mi.beta, pch="", names=paste("col", 1:10), xlab="Column numbers", ylab="Estimated coefficient")
boxplot(res.mi$mi.epsilon, pch="")

library(xtable)
xtable(t(as.matrix(res.lori$alpha*(abs(res.lori$alpha)>1e-2))))
xtable(t(as.matrix(res.lori$beta*(abs(res.lori$beta)>1e-2))))
xtable(t(as.matrix(res.lori$epsilon*(abs(res.lori$epsilon)>1e-2))))

library(plotly)
p <- plot_ly(res.cv$errors, x = ~lambda1, y = ~lambda2, z = ~errors,
             marker=list(size=1)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'lambda1'),
                      yaxis = list(title = 'lambda2'),
                      zaxis = list(title = 'error')))
p


u <- do.call(rbind, lapply(1:m2, function(j) diag(m1)))
u <- cbind(u, sapply(1:m2, function(j) {
  t <- matrix(0,m1,m2)
  t[,j] <- rep(1,m1)
  return(t)
}))
u <- cbind(u,U)

scale(u)

