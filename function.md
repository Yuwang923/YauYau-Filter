# Filter function

## 说明

依赖 `RcppAramdillo` 

## Kalman Filter

### kalmanFiltering

```c++
// [[Rcpp::export]]
mat kalmanFiltering(mat Yobs, vec z0, mat A, mat C, double sigmaQ, double sigmaR) {
    mat Sigma, Sigma_f, S, Kgain, Zmat;
    vec mu, mu_f;
    size_t tIndex, tTotal, p;
    p = z0.n_elem;
    tTotal = Yobs.n_cols;
    mu = z0;
    
    Zmat = zeros<mat>(p, tTotal);
    Sigma = zeros<mat>(p, p);
    
    for(tIndex = 0; tIndex < tTotal; tIndex++){
        // prediction step
        mu_f = A * mu;
        Sigma_f = A * Sigma * A.t();
        Sigma_f.diag() += sigmaQ*sigmaQ;
        
        // compute Kalman gain
        S = C*Sigma_f*C.t();
        S.diag() += sigmaR*sigmaR;
        Kgain = solve(S, C*Sigma_f);
        Kgain = Kgain.t();
        
        // measurement step
        mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f);
        Sigma = Sigma_f - Kgain * C * Sigma_f;
        Zmat.col(tIndex) = mu;
    }
    
    return Zmat;
}
```

### kalmanSmoothing

```c++
// [[Rcpp::export]]
mat kalmanSmoothing(mat Yobs, vec z0, mat A, mat C, double sigmaQ, double sigmaR) {
  // Forward pass: Kalman filtering
  mat Sigma, Sigma_f, S, Kgain, Zmat, mu_smooth, Sigma_smooth;
  vec mu, mu_f;
  size_t tIndex, tTotal, p;
  p = z0.n_elem;
  tTotal = Yobs.n_cols;
  mu = z0;
  
  Zmat = zeros<mat>(p, tTotal);
  Sigma = eye<mat>(p, p) * 1e-6; // Initial covariance set to small value to avoid singular matrix
  std::vector<mat> Sigma_filtered(tTotal);
  std::vector<vec> mu_filtered(tTotal);
  
  for (tIndex = 0; tIndex < tTotal; tIndex++) {
    // Prediction step
    mu_f = A * mu;
    Sigma_f = A * Sigma * A.t();
    Sigma_f.diag() += sigmaQ * sigmaQ;
    
    // Compute Kalman gain
    S = C * Sigma_f * C.t();
    S.diag() += sigmaR * sigmaR;
    Kgain = solve(S, C * Sigma_f);
    Kgain = Kgain.t();
    
    // Measurement update step
    mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f);
    Sigma = Sigma_f - Kgain * C * Sigma_f;
    Zmat.col(tIndex) = mu;
    
    // Store filtered estimates
    mu_filtered[tIndex] = mu;
    Sigma_filtered[tIndex] = Sigma;
  }
  
  // Backward pass: Rauch-Tung-Striebel smoother
  mu_smooth = Zmat;
  Sigma_smooth = Sigma;
  for (tIndex = tTotal - 2; tIndex != (size_t)(-1); tIndex--) {
    mat J = Sigma_filtered[tIndex] * A.t() * inv(Sigma_f);
    
    mu_smooth.col(tIndex) = mu_filtered[tIndex] + J * (mu_smooth.col(tIndex + 1) - A * mu_filtered[tIndex]);
    Sigma_smooth = Sigma_filtered[tIndex] + J * (Sigma_smooth - Sigma_f) * J.t();
    Zmat.col(tIndex) = mu_smooth.col(tIndex);
  }
  
  return Zmat;
}
```

### Rresult

```R

generateData = function(n, A,C, z0, sigmaQ, sigmaR){
    z = z0
    # 4-by-n Actual hidden state
    ZActual = matrix(0, nrow = 4, ncol = n)
    # 2-by-n observation matrix
    Yobs = matrix(0, nrow = 2, ncol = n) 
    # simulate n observations
    for(i in 1:n){
        z = A%*%z + rnorm(4, sd = sigmaQ)
        Yobs[,i] = C %*% z + rnorm(2, sd = sigmaR)
        ZActual[,i] = z
    }
    # return the observation and latent state as a list
    dataList = list(Yobs = Yobs, 
                    ZActual = ZActual)
    
    return(dataList)
}

drawPath = function(Yobs, ZActual, Zhat1, Zhat2){
    plot(t(Yobs),pch = 20, col = "grey", xlab = "", ylab = "")
    lines(t(ZActual[1:2,]), col = "black", lwd = 2)
    lines(t(Zhat1[1:2,]), col = "blue", lwd = 2)
    lines(t(Zhat2[1:2,]), col = "red", lwd = 2)
}

## Parameters
sigmaQ = 0.2 
sigmaR = 0.8
n = 100 # Sample size
delta = 0.5 # step size
A = diag(1, nrow = 4, ncol = 4) # transition matrix
A[1,3] = A[2,4] = delta 
C = diag(1, nrow = 2, ncol = 4) # obs matrix
x0 = c(0,0) # Initial Poisition
v0 = c(1,1) # Initial velocity
z0 = c(x0, v0) # Initial latent state

#### Remove the comments when you finish the function kalmanSmoothing(.)
Rcpp::sourceCpp("homework2/homework2.3/kalmanFiltering.cpp")
Rcpp::sourceCpp("homework2/homework2.3/kalmanSmoothing.cpp")
kdata = generateData(n, A, C, z0, sigmaQ, sigmaR)
Zhat1 = kalmanFiltering(kdata$Yobs, z0, A, C,  sigmaQ, sigmaR)
Zhat2 = kalmanSmoothing(kdata$Yobs, z0, A, C,  sigmaQ, sigmaR)

Yobs = kdata$Yobs
ZActual = kdata$ZActual

jpeg(filename = "~/homework2/homework2.3/Kalman_filtet_smoothing.png", width = 600, height = 600)

par(mar = c(2,2,2,1), mfrow = c(2,1))
plot(t(Yobs),pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Filter")
lines(t(ZActual[1:2,]), col = "grey", lwd = 2)
lines(t(Zhat1[1:2,]), col = "blue", lwd = 2)
plot(t(Yobs),pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Smoothing")
lines(t(ZActual[1:2,]), col = "grey", lwd = 2)
lines(t(Zhat2[1:2,]), col = "blue", lwd = 2)

dev.off()

par(mar = c(2,2,2,1), mfrow = c(2,1))
plot(t(Yobs),pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Filter")
lines(t(ZActual[1:2,]), col = "grey", lwd = 2)
lines(t(Zhat1[1:2,]), col = "blue", lwd = 2)
plot(t(Yobs),pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Smoothing")
lines(t(ZActual[1:2,]), col = "grey", lwd = 2)
lines(t(Zhat2[1:2,]), col = "blue", lwd = 2)

```



## YauYau Filter

### 状态-观测模拟

```c++
\\ 信号-观测模型的模拟

\\ [[Rcpp::export]]
Simulate_State_Obser(){
\\

}
```

```c++
\\ 状态过程的模拟可视化

\\ [[Rcpp::export]]
plot_State(){
\\

}
```

```c++
\\ 观测过程的模拟可视化

\\ [[Rcpp::export]]
plot_Obser(){
\\

}
```

### PDE 离散

### PDE 求解

### 状态估计

### 结果对比
