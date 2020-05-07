####################Real-Data##############################################
# FUNCTION
Adjust_blurring <- function(K_not_normalised, delta){
  K <- K_not_normalised^(delta)
  for (row in c(1:nrow(K))) {
    K[row,] <- K[row,]/sum(K[row,])
  }
  return(K)
}

indexTo2D <- function(L,m){
  # L: index location.
  # m: number of pixels per column.
  i <- ceiling(L/m)
  j <- L%%m
  if(j==0){j=m}
  return(c(i,j))
}

CoordsToIndex <- function(i,j,m){
  # i: x-coordinate.
  # j: y-coordinate.
  # m: number of pixels per column.
  return((i-1)*m+j)
}
Gaussian <- function(delta,im){
  # Input
  # delta: non negative integer representing blurring of the image, higher delta means larger blurring.
  # im: a cimg type.
  K <- matrix(rep(0, dim(im)[1]*dim(im)[2]*dim(im)[1]*dim(im)[2]), nrow = dim(im)[1]*dim(im)[2], ncol = dim(im)[1]*dim(im)[2])
  m <- dim(im)[1]
  for(i in 1:nrow(K)){
    for(j in 1:ncol(K)){
      yx <- ceiling(i/m)
      yy <- i%%m
      if(yy==0){yy=m}
      fx <- ceiling(j/m)
      fy <- j%%m
      if(fy==0){fy=m}
      K[i,j] <- dnorm(x=(abs(yy-fy)+abs(yx-fx)), mean = 0, sd=delta)
    }
    K[i,] <- K[i,]/sum(K[i,])
  }
  return(K)
}
generateFirstOrderL <- function(parameter,height){
  #parameter: number of parameters
  #height: height of images
  L <- matrix(ncol = parameter)
  
  width <- parameter/height
  #width of images
  
  for(i in 1:parameter){
    
    coord <- indexTo2D(i,height)
    #2D-coordinate of i^th parameter
    x <- coord[1]
    y <- coord[2]
    
    
    row <- rep(0,parameter)
    row[i] <- 1
    # The pixel directly below our current pixel. If it exist then assign a value of -1
    if(y+1<=height){
      row[CoordsToIndex(x,y+1,height)] <- -1
      L <- rbind(L,row)
    }
    
    row <- rep(0,parameter)
    row[i] <- 1
    if(x+1<=width){
      row[CoordsToIndex(x+1,y,height)] <- -1
      L <- rbind(L,row)
    }
    
    row <- rep(0,parameter)
    row[i] <- 1
    if(y-1>0 & y-1<=height){
      row[CoordsToIndex(x,y-1,height)] <- -1
      L <- rbind(L,row)
    }
    
    row <- rep(0,parameter)
    row[i] <- 1
    if(x-1>0 & x-1<=width){
      row[CoordsToIndex(x-1,y,height)] <- -1
      L <- rbind(L,row)
    }
  }
  L <- L[rowSums(is.na(L)) != ncol(L), ]
  return(L)
}

# Creating Matrix and Data

delta <- 0.2
K <- Gaussian(delta, Y)
f <- c(Y)
max_f <- max(f)
f <- f/max_f
height <- dim(Y)[2]
width <- dim(Y)[1]
parameter <- length(f)
y_data <- f
sigma <- 0.003763252
####################### Plotting Closed Estimates ####################

data <- as.cimg(matrix(rep(1,parameter) - f, ncol=height))
plot(data)

# Ridge regression with L=I_p
tau <- 0.000001
lambda <- sigma^2/tau^2
estimates_ordinary <- solve(t(K)%*%K + lambda*t(L1)%*%L1,t(K)%*%y_data)
MSE <- mean((estimates_ordinary-f)^2)
par(mfrow=c(1,2))
plot(data, main = paste("Blurred with noise: ", sigma))
estimate2D_or <- matrix(rep(1,length(f))-estimates_ordinary,ncol=height)
estimate2D_or <- as.cimg(estimate2D_or)
plot(estimate2D_or, main = 'Ordinary Ridge Regression')

# Generalised Ridge Regression

tau <- 0.1
lambda <- sigma^2/tau^2
estimates_fused <- solve(t(K)%*%K + lambda*t(L)%*%L,t(K)%*%y_data)
MSE_fuse <- mean((estimates_fused-f)^2)
estimate2D <- matrix(rep(1,length(f))-estimates_fused,ncol=height)
estimate2D <- as.cimg(estimate2D)
par(mfrow=c(1,2))
plot(data, main = paste("Blurred with noise: ", sigma))
plot(estimate2D, main = paste("Generalised Ridsge Regression"))

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
      
      if(log(U) <= new_posterior - old_posterior){
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
      if(log(U) <= tau_new_posterior - tau_old_posterior){
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
RWMCMC_sample_hierarchical <- RWMCMC_hierarchical(y_data,K,sigma, lambda=1,p=parameter,N=2000,g_f=0.02, g_tau = 0.01,q=1, starting_point_tau = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time 
time.taken

burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))

par(mfrow=c(1,2))
plot(as.cimg(matrix(estimate2D, ncol=height)), main = paste("Generalised Ridge Estimator"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(rep(1,parameter)-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),]), ncol = height)), main = 'RWMCMC Hierarchial posterior mean')

par(mfrow=c(1,2))
plot(as.cimg(matrix(1-f, ncol=height)), main = paste("data"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(rep(1,parameter)-colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),]), ncol = height)), main = 'RWMCMC Hierarchial posterior mean')


par(mfrow=c(1,2))
acf(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),1])
acf(RWMCMC_sample_hierarchical$tau[c(burn_in:length(RWMCMC_sample_hierarchical$tau))])
# ERROR



par(mfrow=c(1,1))
error_estimate <- abs(c(estimate2D) -(rep(1,parameter)- colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),])))
error_estimate <- matrix(error_estimate, ncol = height)
plot(as.cimg(error_estimate), main='Difference')

par(mfrow=c(1,1))
error_estimate <- ((1-f) -(colMeans(RWMCMC_sample_hierarchical$sample[c(burn_in:nrow(RWMCMC_sample_hierarchical$sample)),])))^2
error_estimate <- matrix(error_estimate, ncol = height)
plot(as.cimg(error_estimate), main='')

# Posterior Mode

start.time <- Sys.time()
RWMCMC_sample_hierarchical <- RWMCMC_hierarchical(y_data,K,sigma, lambda=1,p=parameter,N=2000,g_f=0.02, g_tau = 0.01,q=1, starting_point_tau = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time 
time.taken

f_hat = RWMCMC_sample_hierarchical$sample[nrow(RWMCMC_sample_hierarchical$sample),]
par(mfrow=c(1,2))
plot(as.cimg(matrix(1-f, ncol=height)), main = paste("Generalised Ridge Estimator"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(1-f_hat, ncol = height)), main = 'RWMCMC Hierarchial posterior mean')

# Gaussian VS Laplace

par(mfrow=c(1,2))
plot(as.cimg(matrix(1-f_Gaus, ncol = height)), main='Gaussian Prior')
plot(as.cimg(matrix(1-f_Lap, ncol = height)), main='Laplacian Prior')

par(mfrow=c(1,1))
plot(as.cimg(matrix((f_Lap-f_Gaus)^2, ncol = height)), main='')

####################### MAP RWMCMC ###########################################

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
    RWMCMC_sample <- matrix(rep(0,p), ncol=p)
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
      
      if(new_posterior >= old_posterior){
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
      plot(RWMCMC_sample[,1], type='l')
    }
    print(c(i/N))
  }
  return(RWMCMC_sample)
}

tau <- 0.008
start.time <- Sys.time()
RWMCMC_sample <- RWMCMC(y_data,K,sigma, tau=tau,p=parameter,N=200,jump_sd=20,q=1)
plot_grid <- function(nx,ny){
  for (x in c(0:nx)) {
    abline(v=x, lty='dotted', col='lightgray')
  }
  for (y in c(0:ny)) {
    abline(h=y, lty = 'dotted', col='lightgray')
  }
}
par(mfrow=c(1,2))
plot(as.cimg(1-Y), main = "Observation")
#plot_grid(dim(Y)[1], dim(Y)[2])
MAP_3 <- RWMCMC_sample[nrow(RWMCMC_sample),]
plot(as.cimg(matrix(1-MAP_3, ncol = height)), main = 'LASSO')
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

par(mfrow=c(1,3))
plot(as.cimg(matrix(1-MAP_1, ncol = height)), main = paste('tau=',1))
plot(as.cimg(matrix(1-MAP_2, ncol = height)), main = paste('tau=',0.015))
plot(as.cimg(matrix(1-MAP_3, ncol = height)), main = paste('tau=',0.008))

MAP_temp <- MAP
MAP_temp[1] <- 1
MAP_temp[MAP_temp<0]=0
par(mfrow=c(1,2))
plot(as.cimg(1-Y), main = "Observation")
plot(as.cimg(matrix(1-MAP_temp, ncol = height)), main = 'LASSO')

####################### Gaussian VS Laplace

f_MAP <- RWMCMC_sample_hierarchical$sample[nrow(RWMCMC_sample_hierarchical$sample),]
tau_MAP <- RWMCMC_sample_hierarchical$tau[length(RWMCMC_sample_hierarchical$tau)]

par(mfrow=c(1,2))
plot(as.cimg(matrix(f, ncol=height)), main = paste("data"))
#plot(im.grayscale.resize, main = "Truth")
burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))
plot(as.cimg(matrix(1-f_MAP, ncol = height)), main = 'RWMCMC Hierarchial posterior mean')

###################### Noise estimation

noise_estimate <- c(Y[c(13:16),35], 
                    Y[c(12:17),36], 
                    Y[c(11:18),37], 
                    Y[c(10:19),38], 
                    Y[c(9:20),39], 
                    Y[c(9:20),40], 
                    Y[c(9:20),41], 
                    Y[c(9:20),42], 
                    Y[c(9:20),43], 
                    Y[c(10:19),44], 
                    Y[c(11:18),45], 
                    Y[c(12:17),46], 
                    Y[c(13:17),47])
noise_estimate <- noise_estimate/max_f
var(noise_estimate)
####################### Hierarchical MAP RWMCMC ##########################

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
      
      if(log(U) <= new_posterior - old_posterior){
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
      if(log(U) <= tau_new_posterior - tau_old_posterior){
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
        gamma_f <- gamma_f*(rgt_f + 0.23)/0.46
        gamma_tau <- gamma_tau*(rgt_tau + 0.23)/0.46
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
RWMCMC_sample_hierarchical <- RWMCMC_hierarchical(y_data,K,sigma, lambda=10^(6),p=parameter,N=10000,g_f=0.1, g_tau = 0.1,q=2, starting_point_tau = 3)
end.time <- Sys.time()
time.taken <- end.time - start.time 
time.taken

par(mfrow=c(1,2))
ACF <- acf(RWMCMC_sample_hierarchical$sample[,1], main='top left corner pixel')
acf(RWMCMC_sample_hierarchical$tau, main='tau')

sum(ACF$acf)

ACF_sum <- 0
for (i in c(1:ncol(RWMCMC_sample_hierarchical$sample))) {
  ACF <- acf(RWMCMC_sample_hierarchical$sample[,i], main='top left corner pixel')
  ACF_sum <- ACF_sum + sum(ACF$acf)
}
ACF_sum/ncol(RWMCMC_sample_hierarchical$sample)

burn_in <- floor(0.3*nrow(RWMCMC_sample_hierarchical$sample))

f_MAP <- RWMCMC_sample_hierarchical$sample[nrow(RWMCMC_sample_hierarchical$sample),]
tau_MAP <- RWMCMC_sample_hierarchical$tau[length(RWMCMC_sample_hierarchical$tau)]
tau_MAP
f_MAP[f_MAP<0]=0
par(mfrow=c(1,1))
plot(as.cimg(1-Y), main = "Observation")
plot(as.cimg(matrix(1-f_MAP, ncol = height)), main = '')
#par(mfrow=c(1,1))
plot(as.cimg(matrix((f-f_MAP)^2, ncol = height)), main = '')
save(f_MAP, file='mouse2_lambda_10_6.RData')