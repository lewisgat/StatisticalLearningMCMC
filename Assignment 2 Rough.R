##############################
### Obtaining initial Data ###
##############################

library(datasets)
library(MASS)
library(mvtnorm)
data(iris)
index_versicolor = which(iris[,5]=='versicolor')
index_virginica = which(iris[,5]=='virginica')
index_setosa = which(iris[,5]=='setosa')
iris[,5] = 1
iris[index_virginica,5] = 2
iris[index_setosa,5] = 3
iris[,5] = factor(iris[,5])
x <- iris[1:150, c("Sepal.Length", "Sepal.Width", "Species")]
load(file = "index_iris_2182973.Rdata")
number_train = 149
x_train = as.matrix(x[index[1:number_train],1:2])
y_train = x[index[1:number_train],3]
x_test = as.matrix(x[index[number_train+1],1:2])
y_test = x[index[number_train+1],3]

X_train = cbind(matrix(1, dim(x_train)[1], 1), x_train)

#####################################
### Creating Likelihood Functions ###
#####################################


sigmoid_func <- function(X_vec, beta_vec){
  
  X_m <- matrix(X_vec, nrow = 1, ncol = 3)
  
  beta_M <- matrix(beta_vec, nrow = 3, ncol = 1)
  
  Z <- X_m %*% beta_M  #defining sigmoid function
  return(Z)
}

log_likelihood <- function(X, y, beta){
  
  l <- 0
  n <- dim(X)[1]
  
  for (i in 1:n){
    
    Z_1 <- sigmoid_func(X[i,], beta[1:3])
    Z_2  <- sigmoid_func(X[i,], beta[4:6])
    Z_3 <- sigmoid_func(X[i,], beta[7:9])
    
    l <- l + (y[i] == 1)*Z_1 + (y[i]==2)*Z_2
            + (y[i] == 3)*Z_3 - log(exp(Z_1) + exp(Z_2) + exp(Z_3))
  }
  return(l) }


##########################
### First MH Algorithm ###
##########################



Metropolis_Hastings <- function(X, y, iterations, beta0, sigma_q) {
  acc <- 0
  betas <- array(NA, dim = c(iterations+1, 9))
  betas[1,] <- beta0
  
  for (i in 1:iterations) {
    print(i)
    beta_star <- mvrnorm(n = 1, betas[i,], sigma_q)
    
    
    p <- (log_likelihood(X, y, beta_star) + 
           log(dmvnorm(beta_star, mean = matrix(0,9,1), sigma = diag(9) * 100)) + 
            log(dmvnorm(betas[i, ], mean = beta_star, sigma = sigma_q)) ) - 
      (log_likelihood(X, y, betas[i, ]) + 
          log(dmvnorm(betas[i, ], mean = matrix(0,9,1), sigma = diag(9) * 100)) + 
           log(dmvnorm(beta_star, mean = betas[i,], sigma = sigma_q)) )
    
    
    if (runif(1) < min(1, p)) {
      betas[i+1,] <- beta_star
      acc <- acc + 1
    } 
    else {
      betas[i+1,] <- betas[i,]
    }
      
  }
  
  
  print(acc/iterations)
  return(betas)  
}

#matrix(0.25, nrow = 9, ncol = 9) +
#sigma_q = (2.38**2) * (diag(9)*0.05) / 9

# r <-  matrix(rnorm(81,0,0.05),nrow = 9, ncol = 9)
# sigma_1 <- t(r) %*% r
# 
# for (i in 1:9){
#   sigma_1[i,i] <- 0.025
# 
# }

sigma_1 <- diag(9)*0.025


beta0 = c(1,1,1,1,1,1,1,1,1)

betas_gen <- Metropolis_Hastings(X_train, y_train, iterations = 200000, beta0, sigma_1)

plot(betas_gen[,5], xlim = c(1,1001))

