###################loading packages########################################
library('imager')
library('MASS')

###################loading image###########################################
file <- "./bottom_left.png"
im <- load.image(file)
plot(im)

# Processing image and rescaling image
im.rgb <- rm.alpha(im)
im.grayscale <- grayscale(im.rgb)
im.grayscale.resize <- resize(im.grayscale, size_x= 40, size_y=30)
plot(im.grayscale.resize)

# Vectorising image
f <- c(im.grayscale.resize) # truth
####################Simulation##############################################
sigma <- 0.02 # Variance of noise
load("./simulated/K.RData") # Blur matrix
load("./simulated/L.RData") # Penalisation matrix
mean <- rep(0,length(f))
Cov <- sigma^2*diag(length(f))
blurred <- K%*%f
noise <- mvrnorm(n=1,mu=mean, Cov)
y <- blurred + noise
y <- 1-y
obs <- matrix(1-y,nrow=dim(im.grayscale.resize)[1])
obs <- as.cimg(obs)
par(mfrow=c(1,2))
plot(im.grayscale.resize, main = "Orginal")
plot(as.cimg(matrix(obs,nrow=dim(im.grayscale.resize)[1])), main = "")

height <- dim(obs)[2]
width <- dim(obs)[1]
parameter <- length(f)
L1 <- diag(parameter)

#####################Ridge regression with L=I_p#############################
tau <- 0.0988687
tau <- 0.06
lambda <- sigma^2/tau^2
estimates_ordinary <- solve(t(K)%*%K + lambda*t(L1)%*%L1,t(K)%*%y)
MSE <- mean((estimates_ordinary-f)^2)
par(mfrow=c(1,3))
plot(obs, main = paste("Blurred with noise: ", sigma))
plot(im.grayscale.resize, main = "Orginal")
estimate2D_or <- matrix(1-estimates_ordinary,ncol=height)
estimate2D_or <- as.cimg(estimate2D_or)
plot(estimate2D_or, main = paste('Ordinary Ridge Regression'))

######################Generalised Ridge Regression###########################

tau <- 0.04774989
tau <- 0.1
lambda <- sigma^2/tau^2
estimates_fused <- solve(t(K)%*%K + lambda*t(L)%*%L,t(K)%*%y)
MSE_fuse <- sum((estimates_fused-f)^2)
estimate2D <- matrix(rep(1,length(f))-estimates_fused,ncol=height)
estimate2D <- as.cimg(estimate2D)
par(mfrow=c(1,3))
plot(obs, main = paste("Blurred with noise: ", sigma))
plot(im.grayscale.resize, main = "Orginal")
plot(estimate2D, main = paste("Generalised Ridsge Regression"))

# comparison between L and I

par(mfrow=c(1,2))
plot(estimate2D_or, main = 'Zero Order')
plot(estimate2D, main = 'First Order')

####################### Optimal tau ##################################

find_MSE <- function(input, truth){
  return(sum((input-truth)^2))
}
j <- 1
for (i in c(1:100)) {
  tau[j] <- i*10^(-2)
  lambda <- sigma^2/tau[j]^2
  estimates_ordinary <- solve(t(K)%*%K + lambda*t(L1)%*%L1,t(K)%*%y)
  MSE_new <- find_MSE(1-estimates_ordinary,f)
  MSE[j] <- MSE_new
  j <- j+1
}
par(mfrow=c(1,1))
plot(x=tau, y=MSE, type='l')
####################### RWMCMC ###############################################

RWMCMC <- function(y,K,sigma,tau,p,N,jump_sd, q,starting_point=NA){
  #INPUT:
  # y (array)       : Array of observation
  # K (matrix)      : Blurring Matrix
  # sigma (constant): Variance of observation
  # tau (constant)  : Prior variance of parameters
  # p (integer)     : Number of parameters
  # N (integer)     : Number of samples
  # jump_sd         : variance of jump distribution
  # starting_point  : a vector of starting point, if known.
  #OUTPUT:
  # RWMCMC (matrix) : A matrix of samples where j th column is the j th parameter
  #                   and the i th row is the i th sample.
  
  if(is.na(starting_point)){
    RWMCMC_sample <- matrix(rep(1,p), ncol=p)
  }else{
    RWMCMC_sample <- matrix(starting_point, ncol=p)  
  }
  
  old_f <- RWMCMC_sample[1,]
  old_Kf <- K%*%old_f
  old_Lf <- L%*%old_f
  old_like <- -0.5*(sigma)^(-2)*sum((y-old_Kf)^2)
  old_posterior <- old_like -0.5*(tau)^(-2)*sum(abs(old_Lf)^q)
  for (i in c(2:N)) {
    for (j in c(1:p)) {
      epsilon <- rnorm(1, 0, sd=jump_sd)
      U <- runif(1,0,1)
      
      new_like <- old_like + (sigma)^(-2)*sum((y-old_Kf)*K[,j])*epsilon -0.5*(sigma)^(-2)*epsilon^2*sum(K[,j]^2)
      
      new_Lf <- old_Lf + L[,j]*epsilon
      
      new_posterior <- new_like -0.5*(tau)^(-2)*sum(abs(new_Lf)^q)
      
      if(log(U) <= new_posterior - old_posterior){
        new_Kf <- old_Kf + K[,j]*epsilon
        
        old_f[j] <- old_f[j] + epsilon 
        old_posterior <- new_posterior
        old_like <- new_like
        old_Kf <- new_Kf
        old_Lf <- new_Lf
      }
    }
    RWMCMC_sample <- rbind(RWMCMC_sample, old_f)
    if(i%%10==0){
      par(mfrow=c(1,1))
      plot(RWMCMC_sample[,1], type='l', ylab='')
    }
    print(c(i/N))
  }
  return(RWMCMC_sample)
}

tau <- 0.15
start.time <- Sys.time()
RWMCMC_sample <- RWMCMC(y,K,sigma, tau=tau,p=parameter,N=1000,jump_sd=0.008,q=1)
par(mfrow=c(1,3))
plot(estimate2D, main = paste("Generalised Ridge Estimator"))
plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample))
plot(as.cimg(matrix(1-RWMCMC_sample[nrow(RWMCMC_sample),], ncol = height)), main = 'LASSO')
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

####################### Hierarchical RWMCMC ##########################

RWMCMC_hierarchical <- function(y,K,sigma,lambda,p,N,g_f, g_tau, q,starting_point=NA, starting_point_tau = NA){
  #INPUT:
  # y (array)       : Array of observation
  # K (matrix)      : Blurring Matrix
  # sigma (constant): Variance of observation
  # tau (constant)  : Prior variance of parameters
  # p (integer)     : Number of parameters
  # N (integer)     : Number of samples
  # jump_sd         : variance of jump distribution
  # starting_point  : a vector of starting point, if known.
  #OUTPUT:
  # RWMCMC (matrix) : A matrix of samples where j th column is the j th parameter
  #                   and the i th row is the i th sample.
  
  if(is.na(starting_point)){
    RWMCMC_sample <- matrix(rep(0,p), ncol=p)
  }else{
    RWMCMC_sample <- matrix(starting_point, ncol=p)  
  }
  if (is.na(starting_point_tau)) {
    tau <- 3
  }else{
    tau <- starting_point_tau
  }
  gamma_f <- g_f
  gamma_tau <- g_tau
  
  accept_f <- 0
  total_f <- 0
  accept_tau <- 0
  total_tau <- 0
  
  burn_in <- floor(0.3*N)
  tau_sample <- c(tau)
  old_f <- RWMCMC_sample[1,]
  old_Kf <- K%*%old_f
  old_Lf <- L%*%old_f
  old_like <- -0.5*(sigma)^(-2)*sum((y-old_Kf)^2)
  old_posterior <- old_like -0.5*(tau)^(-2)*sum(abs(old_Lf)^q)  - parameter*log(2*pi*tau^2) - lambda*tau^2
  tau_old_posterior <- old_like + -0.5*(tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
  for (i in c(2:N)) {
    for (j in c(1:p)) {
      epsilon <- rnorm(1, 0, sd=gamma_f)
      U <- runif(1,0,1)
      
      new_like <- old_like + (sigma)^(-2)*sum((y-old_Kf)*K[,j])*epsilon -0.5*(sigma)^(-2)*epsilon^2*sum(K[,j]^2)
      
      new_Lf <- old_Lf + L[,j]*epsilon
      
      new_posterior <- new_like -0.5*(tau)^(-2)*sum(abs(new_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
      
      if(new_posterior >= old_posterior){
        new_Kf <- old_Kf + K[,j]*epsilon
        
        old_f[j] <- old_f[j] + epsilon 
        old_posterior <- new_posterior
        old_like <- new_like
        old_Kf <- new_Kf
        old_Lf <- new_Lf
        accept_f <- accept_f + 1
      }
      total_f <- total_f + 1
    }
    
    new_tau <- rnorm(1, tau, sd=gamma_tau)
    if(new_tau > 0){
      tau_old_posterior <- old_like -0.5*(tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
      tau_new_posterior <- old_like -0.5*(new_tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*new_tau^2) - lambda*new_tau^2
      U <- runif(1,0,1)
      if(tau_new_posterior >= tau_old_posterior){
        tau <- new_tau
        accept_tau <- accept_tau + 1
        #print(paste('proposed tau: ', new_tau,', accepted'))
      }
    }
    tau_sample <- c(tau_sample,tau)
    total_tau <- total_tau + 1
    RWMCMC_sample <- rbind(RWMCMC_sample, old_f)
    if(i%%10==0){
      par(mfrow=c(2,1))
      plot(RWMCMC_sample[,1], type='l')
      plot(tau_sample, type='l')
      if(i<=burn_in && i%%floor(burn_in/6)==0){
        rgt_f <- accept_f/total_f
        rgt_tau <- accept_tau/total_tau
        gamma_f <- gamma_f*(rgt_f + 0.4)/0.8
        gamma_tau <- gamma_tau*(rgt_tau + 0.4)/0.8
        print(paste(' f variance: ', gamma_f, ' tau variance: ', gamma_tau))
        accept_tau <- 0
        total_tau <- 0
        accept_f <- 0
        accept_f <- 0
      }
    }
    print(c(i/N))
  }
  return(list('sample'=RWMCMC_sample,'tau'=tau_sample))
}

start.time <- Sys.time()
RWMCMC_sample_hierarchical <- RWMCMC_hierarchical(y,K,sigma, lambda=10^3,p=parameter,N=300,g_f=0.05, g_tau = 0.1,q=1, starting_point_tau = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time 
time.taken

burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))

mean(RWMCMC_sample_hierarchical$tau[c(burn_in:length(RWMCMC_sample_hierarchical$tau))])

par(mfrow=c(1,2))
plot(as.cimg(matrix(estimate2D, ncol=height)), main = '')
#plot(im.grayscale.resize, main='')
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
estimator <- colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),])
plot(as.cimg(matrix(rep(1,parameter)-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),]), ncol = height)), main = '')

par(mfrow=c(1,2))
plot(as.cimg(matrix(f, ncol=height)), main = paste("data"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(rep(1,parameter)-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),]), ncol = height)), main = 'RWMCMC Hierarchial posterior mean')


par(mfrow=c(1,2))
acf(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),1], main='Top Left Corner Pixel')
acf(RWMCMC_sample_hierarchical$tau[c(burn_in:length(RWMCMC_sample_hierarchical$tau))], main=paste('tau'))


tcp <- acf(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),1], main='Top Left Corner Pixel')
print(sum(abs(tcp$acf)))
total <- 0
for (col in c(1:ncol(RWMCMC_sample_hierarchical$sample))) {
  temp <- acf(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),col], main='Top Left Corner Pixel')
  total <- total + sum(abs(temp$acf))
}
print(total / 1200)

# Gauss ERROR
max <- 1
par(mfrow=c(1,1))
error_estimate <- (f -(1-estimates_fused))^2
print(sum(error_estimate))
max <- max(error_estimate)
error_estimate[1] <- max
error_estimate <- matrix(error_estimate, ncol = height)
plot(as.cimg(error_estimate), main='Gaussian Error')

# Laplace ERROR
par(mfrow=c(1,1))
error_estimate <- (f -(1-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),])))^2
print(sum(error_estimate))
error_estimate[1] <- max
error_estimate <- matrix(error_estimate, ncol = height)
plot(as.cimg(error_estimate), main='Laplace Error')

####################### Gaussian VS Laplace

f_MAP <- RWMCMC_sample_hierarchical$sample[nrow(RWMCMC_sample_hierarchical$sample),]
tau_MAP <- RWMCMC_sample_hierarchical$tau[length(RWMCMC_sample_hierarchical$tau)]
f_MAP[f_MAP<0]=0
par(mfrow=c(1,1))
plot(as.cimg(matrix(1-Laplace, ncol=height)), main = paste(""))
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(1-f_MAP, ncol = height)), main = '')
plot(im.grayscale.resize, main = "")

sum((f-(1-f_MAP))^2)

####################### Hierarchical RWMCMC 3 with prior on Blurring##########################

Adjust_blurring <- function(K_not_normalised, delta){
  K <- K_not_normalised^(delta)
  for (row in c(1:nrow(K))) {
    K[row,] <- K[row,]/sum(K[row,])
  }
  return(K)
}

RWMCMC_hierarchical <- function(y,K_not_normalised,sigma,lambda,p,N,jump_sd, q,starting_point=NA){
  #INPUT:
  # y (array)       : Array of observation
  # K (matrix)      : Blurring Matrix
  # sigma (constant): Variance of observation
  # tau (constant)  : Prior variance of parameters
  # p (integer)     : Number of parameters
  # N (integer)     : Number of samples
  # jump_sd         : variance of jump distribution
  # starting_point  : a vector of starting point, if known.
  #OUTPUT:
  # RWMCMC (matrix) : A matrix of samples where j th column is the j th parameter
  #                   and the i th row is the i th sample.
  
  if(is.na(starting_point)){
    RWMCMC_sample <- matrix(rep(0,p), ncol=p)
  }else{
    RWMCMC_sample <- matrix(starting_point, ncol=p)  
  }
  
  
  gamma_f <- 0.05
  gamma_tau <- 0.005
  gamma_delta <- 0.002
  
  accept_f <- 0
  total_f <- 0
  accept_tau <- 0
  total_tau <- 0
  accept_delta <- 0
  total_delta <- 0
  
  burn_in <- floor(0.3*N)
  tau <- 1
  tau_sample <- c(tau)
  
  delta <- 1
  K <- Adjust_blurring(K_not_normalised, 1/delta)
  delta_sample <- c(delta)
  
  old_f <- RWMCMC_sample[1,]
  old_Kf <- K%*%old_f
  old_Lf <- L%*%old_f
  old_like <- -0.5*(sigma)^(-2)*sum((y-old_Kf)^2)
  old_posterior <- old_like -0.5*(tau)^(-2)*sum(abs(old_Lf)^q)  - parameter*log(2*pi*tau^2) - lambda*tau^2
  tau_old_posterior <- old_like + -0.5*(tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
  for (i in c(2:N)) {
    for (j in c(1:p)) {
      epsilon <- rnorm(1, 0, sd=gamma_f)
      U <- runif(1,0,1)
      
      new_like <- old_like + (sigma)^(-2)*sum((y-old_Kf)*K[,j])*epsilon -0.5*(sigma)^(-2)*epsilon^2*sum(K[,j]^2)
      
      new_Lf <- old_Lf + L[,j]*epsilon
      
      new_posterior <- new_like -0.5*(tau)^(-2)*sum(abs(new_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
      if(log(U)<= new_posterior - old_posterior){
        new_Kf <- old_Kf + K[,j]*epsilon
        
        old_f[j] <- old_f[j] + epsilon 
        old_posterior <- new_posterior
        old_like <- new_like
        old_Kf <- new_Kf
        old_Lf <- new_Lf
        accept_f <- accept_f + 1
      }
      total_f <- total_f + 1
    }
    
    new_tau <- rnorm(1, tau, sd=gamma_tau)
    if(new_tau > 0){
      tau_old_posterior <- old_like -0.5*(tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*tau^2) - lambda*tau^2
      tau_new_posterior <- old_like -0.5*(new_tau)^(-2)*sum(abs(old_Lf)^q) - parameter*log(2*pi*new_tau^2) - lambda*new_tau^2
      U <- runif(1,0,1)
      if(log(U)<= tau_new_posterior - tau_old_posterior){
        tau <- new_tau
        accept_tau <- accept_tau + 1
        #print(paste('proposed tau: ', new_tau,', accepted'))
      }
    }
    tau_sample <- c(tau_sample,tau)
    total_tau <- total_tau + 1
    
    new_delta <- rnorm(1, delta, gamma_delta)
    
    if(new_delta>0){
      
      new_K <- Adjust_blurring(K_not_normalised, delta/new_delta)
      
      delta_old_posterior <- old_like - 100*delta
      
      new_like <- -0.5*(sigma)^(-2)*sum((y-new_K%*%old_f)^2)
      
      delta_new_posterior <- new_like - 100*new_delta
      
      #print(paste('new delta: ',new_delta,', old delta: ', delta,', new p: ',delta_new_posterior,', old p:' ,delta_old_posterior))
      
      U <- runif(1,0,1)
      
      if(log(U)<= delta_new_posterior - delta_old_posterior){
        K <- new_K
        
        delta <- new_delta
        
        accept_delta <- accept_delta + 1
        
        old_like <- new_like
      }
    }
    total_delta <- total_delta + 1
    delta_sample <- c(delta_sample, delta)
    
    RWMCMC_sample <- rbind(RWMCMC_sample, old_f)
    if(i%%10==0){
      par(mfrow=c(3,1))
      plot(RWMCMC_sample[,1], type='l')
      plot(tau_sample, type='l')
      plot(delta_sample, type='l')
      # if(i<=burn_in && i%%floor(burn_in/6)==0){
      #   rgt_f <- accept_f/total_f
      #   rgt_tau <- accept_tau/total_tau
      #   rgt_delta <- accept_delta/total_delta
      #   gamma_f <- gamma_f*(rgt_f + 0.234)/0.468
      #   gamma_tau <- gamma_tau*(rgt_tau + 0.234)/0.468
      #   gamma_delta <- gamma_delta*(rgt_delta + 0.234)/0.468
      #   print(paste(' f variance: ', gamma_f, ' tau variance: ', gamma_tau))
      #   accept_tau <- 0
      #   accept_f <- 0
      #   accept_delta <- 0
      #   total_tau <- 0
      #   total_f <- 0
      #   total_delta <- 0
      # }
    }
    print(c(i/N))
  }
  return(list('sample'=RWMCMC_sample,'tau'=tau_sample, 'delta'=delta_sample))
}

start.time <- Sys.time()
RWMCMC_sample_hierarchical <- RWMCMC_hierarchical(y,K_not_normalised,sigma, lambda=10,p=parameter,N=1000,jump_sd=0.05,q=1)
par(mfrow=c(1,2))
plot(estimate2D, main = paste("Generalised Ridge Estimator"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(rep(1,parameter)-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),]), ncol = height)), main = 'RWMCMC Hierarchial posterior mean')
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
par(mfrow=c(1,2))
acf(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),1])
acf(RWMCMC_sample_hierarchical$tau[c(burn_in:length(RWMCMC_sample_hierarchical$tau))])
# ERROR
par(mfrow=c(1,1))
error_estimate <- (f -(rep(1,parameter)- colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),])))^2
error_estimate <- matrix(error_estimate, ncol = height)
plot(as.cimg(error_estimate), main='error')


