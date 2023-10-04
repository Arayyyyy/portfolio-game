library(CVXR)

GMVP_portfolio_fun <- function(dataset, ...) {
  X <- diff(log(dataset$adjusted))[-1]  # compute log returns
  Sigma <- cov(X)  # compute SCM
  # design GMVP
  w <- solve(Sigma, rep(1, nrow(Sigma)))
  w <- abs(w)/sum(abs(w))
  return(w)
}

VaR<-function(dataset,...){
    
    X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
    mu <- colMeans(X)  # compute mean vector
    Sigma <- cov(X)  # compute the SCM
    # design mean-variance portfolio
    N<-nrow(Sigma)
    w <- rep(1/N,N)
    
    n<-2000
    VaR_list<-c()
    for (i in 1:n){
      w<-c(runif(N,0,1))
      w<-w/sum(w)
      VaR<-quantile(diff(log(dataset$adjusted%*%w))[-1],0.01)
      
      VaR_list<-append(VaR_list,VaR)
      
      if (VaR == min(VaR_list))
      {
        w_max=w
      }
    }
    
    return(w_max)
    
}

simple<-function(dataset,...){
  
  X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  N<-nrow(Sigma)
  w<-c()
  for(i in 1:N){
    if(mu[i]>0)
    {
      w[i]=mu[i]
    }else
      if(mu[i]<0)
      {
        w[i]=0
      }
  }
  w<-w/sum(w)
  return(w)
}

simple_pro<-function(dataset,...)
{
  X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  N<-nrow(Sigma)
  w<-c()
  for(i in 1:N){
    w[i]=mu[i]/Sigma[i,i]
  }
  
  
  w<-w/sum(abs(w))
  return(w)
}

min_max<-function(dataset,...){
  X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  N<-nrow(Sigma)
  w<-c()
  Max = max(mu)
  Min = min(mu)
  
  for(i in 1:N){
    w[i]=(mu[i]-Min)/(Max-Min)
  }
  w<-w/sum(w)
  return(w)
}

MDCP <- function(dataset,...) {
  X <- as.matrix(diff(log(dataset$adjusted))[-1])
  Sigma<-cov(X)
  N<-nrow(Sigma)
  C <- diag(1/sqrt(diag(Sigma))) %*% Sigma %*% diag(1/sqrt(diag(Sigma)))
  colnames(C) <- colnames(Sigma)
  n<-5000
  pro_list<-c()
  for (i in 1:n){
    w<-c(runif(N,0,1))
    w<-w/sum(w)
    pro<-t(w)%*%C%*%w
    
    pro_list <-append(pro_list,pro)
    if (pro == min(pro_list))
    {
      w_max=w
    }
  }
  return(w_max)
}


MDP <-function(dataset,...){
  
  X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  N<-nrow(Sigma)
  w <- Variable(nrow(Sigma))
  Max = max(mu)
  Min = min(mu)
  
  prob <- Problem(Minimize(quad_form(w, Sigma)),
                  constraints = list( t( sqrt(diag(Sigma))) %*% w == 1))
  result <- CVXR::solve(prob)
  w <- as.vector(result$getValue(w)/sum(abs(result$getValue(w))))
  return(w)
}


rho_LW <- function(X, Sigma_scm, Sigma_T) {
  X <- scale(X, center = TRUE, scale = FALSE)
  T <- nrow(X)
  # tmp <- 0
  # for (t in 1:T)
  #  tmp <- tmp + norm(Sigma_scm - X[t, ] %o% X[t, ], "F")^2
  tmp <- (2-T)*norm(Sigma_scm, "F")^2 + sum(apply(X, 1, norm, "2")^4)
  rho <- min(1, (1/T^2) * tmp / norm(Sigma_scm - Sigma_T, "F")^2)
  return(rho)
}

rho_maxSR <- function(mu_sm, Sigma_scm, Sigma_T, T) {
  W <- diag(T) - (1/T)*matrix(1, T, T)
  rho1_sweep <- exp(seq(-10, 10, length.out = 100))
  obj <- NULL
  for (rho1 in rho1_sweep) {
    Sigma_sh <- rho1*Sigma_T + Sigma_scm
    D <- (1/T)* sum(diag(Sigma_scm %*% solve(Sigma_sh)))
    delta <- D/(1-D)
    b <- T/sum(diag(W %*% solve((diag(T)+delta*W) %*% (diag(T)+delta*W))))
    inv_S_sh_mu <- solve(Sigma_sh, mu_sm)
    num <- mu_sm %*% inv_S_sh_mu - delta
    den <- sqrt(b * inv_S_sh_mu %*% Sigma_scm %*% inv_S_sh_mu)
    obj <- c(obj, num/den)
  }
  #plot(obj)
  i_max <- which.max(obj)
  rho1 <- rho1_sweep[i_max]
  return(rho1/(1 + rho1))
}

GMVP_shrinkage <-function(dataset,...){
  
  X <- as.matrix(diff(log(dataset$adjusted))[-1])  # compute log returns
  mu <- colMeans(X)  # compute mean vector
  Sigma <- cov(X)  # compute the SCM
  # design mean-variance portfolio
  N<-nrow(Sigma)
  w <- Variable(nrow(Sigma))
  T <-nrow(X)
  T_I <- sum(diag(Sigma))/N * diag(N)
  T_D <- diag(diag(Sigma))
  
  rho <- rho_LW(X, Sigma, T_I)
  Sigma_LW_I <- (1-rho)*Sigma + rho*T_I
  rho <- rho_LW(X, Sigma, T_D)
  Sigma_LW_D <- (1-rho)*Sigma + rho*T_D
  
  
  w <- solve(Sigma_LW_D, rep(1, nrow(Sigma_LW_D)))
  w <- (1.5*w)/sum(abs(w))
  return(w)
  
  
  
}

