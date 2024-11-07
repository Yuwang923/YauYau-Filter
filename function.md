# Filter function

## 说明

依赖 `RcppAramdillo` 

## Kalman Filter

### kalmanFiltering

```c++
// [[Rcpp::export]]
mat kalmanFiltering(mat Yobs, vec z0, mat A, mat C, double sigmaQ, double sigmaR) {
    // 定义函数 kalmanFiltering，用于实现卡尔曼滤波器
    // 输入参数：
    // Yobs - 观测数据矩阵，每列代表一个时间点的观测值
    // z0 - 初始状态向量
    // A - 状态转移矩阵，用于预测下一个状态
    // C - 观测矩阵，用于将状态映射到观测空间
    // sigmaQ - 过程噪声标准差
    // sigmaR - 观测噪声标准差
    
    mat Sigma, Sigma_f, S, Kgain, Zmat;
    vec mu, mu_f;
    size_t tIndex, tTotal, p;
    
    // 获取状态向量的维度 p
    p = z0.n_elem;
    // 获取观测数据的时间步数 tTotal
    tTotal = Yobs.n_cols;
    
    // 初始化 mu 为初始状态向量 z0
    mu = z0;
    
    // 初始化存储状态估计值的矩阵 Zmat，大小为 p x tTotal
    Zmat = zeros<mat>(p, tTotal);
    // 初始化状态协方差矩阵 Sigma 为零矩阵，大小为 p x p
    Sigma = zeros<mat>(p, p);
    
    // 循环遍历每个时间步 tIndex
    for(tIndex = 0; tIndex < tTotal; tIndex++){
        // 预测步骤
        // 预测状态均值 mu_f = A * mu
        mu_f = A * mu;
        // 预测状态协方差 Sigma_f = A * Sigma * A^T + Q，其中 Q 表示过程噪声协方差
        Sigma_f = A * Sigma * A.t();
        Sigma_f.diag() += sigmaQ * sigmaQ;
        
        // 计算卡尔曼增益
        // 计算观测预测协方差 S = C * Sigma_f * C^T + R，其中 R 表示观测噪声协方差
        S = C * Sigma_f * C.t();
        S.diag() += sigmaR * sigmaR;
        // 计算卡尔曼增益 Kgain = Sigma_f * C^T * S^{-1}
        Kgain = solve(S, C * Sigma_f);
        Kgain = Kgain.t();
        
        // 测量更新步骤
        // 更新状态均值 mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f)
        mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f);
        // 更新状态协方差 Sigma = Sigma_f - Kgain * C * Sigma_f
        Sigma = Sigma_f - Kgain * C * Sigma_f;
        
        // 存储当前时间步的状态估计值到 Zmat 中
        Zmat.col(tIndex) = mu;
    }
    
    // 返回状态估计矩阵 Zmat
    return Zmat;
}
```

### kalmanSmoothing

```c++
// [[Rcpp::export]]
mat kalmanSmoothing(mat Yobs, vec z0, mat A, mat C, double sigmaQ, double sigmaR) {
  // 定义函数 kalmanSmoothing，用于实现卡尔曼平滑器
  // 输入参数：
  // Yobs - 观测数据矩阵，每列代表一个时间点的观测值
  // z0 - 初始状态向量
  // A - 状态转移矩阵，用于预测下一个状态
  // C - 观测矩阵，用于将状态映射到观测空间
  // sigmaQ - 过程噪声标准差
  // sigmaR - 观测噪声标准差

  // 前向过程：卡尔曼滤波
  mat Sigma, Sigma_f, S, Kgain, Zmat, mu_smooth, Sigma_smooth;
  vec mu, mu_f;
  size_t tIndex, tTotal, p;

  // 获取状态向量的维度 p
  p = z0.n_elem;
  // 获取观测数据的时间步数 tTotal
  tTotal = Yobs.n_cols;

  // 初始化 mu 为初始状态向量 z0
  mu = z0;

  // 初始化存储状态估计值的矩阵 Zmat，大小为 p x tTotal
  Zmat = zeros<mat>(p, tTotal);
  // 初始化状态协方差矩阵 Sigma，大小为 p x p，初值设为很小的值，避免奇异矩阵
  Sigma = eye<mat>(p, p) * 1e-6;

  // 创建用于存储滤波结果的向量
  std::vector<mat> Sigma_filtered(tTotal);
  std::vector<vec> mu_filtered(tTotal);

  // 循环遍历每个时间步 tIndex
  for (tIndex = 0; tIndex < tTotal; tIndex++) {
    // 预测步骤
    // 预测状态均值 mu_f = A * mu
    mu_f = A * mu;
    // 预测状态协方差 Sigma_f = A * Sigma * A^T + Q，其中 Q 表示过程噪声协方差
    Sigma_f = A * Sigma * A.t();
    Sigma_f.diag() += sigmaQ * sigmaQ;

    // 计算卡尔曼增益
    // 计算观测预测协方差 S = C * Sigma_f * C^T + R，其中 R 表示观测噪声协方差
    S = C * Sigma_f * C.t();
    S.diag() += sigmaR * sigmaR;
    // 计算卡尔曼增益 Kgain = Sigma_f * C^T * S^{-1}
    Kgain = solve(S, C * Sigma_f);
    Kgain = Kgain.t();

    // 测量更新步骤
    // 更新状态均值 mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f)
    mu = mu_f + Kgain * (Yobs.col(tIndex) - C * mu_f);
    // 更新状态协方差 Sigma = Sigma_f - Kgain * C * Sigma_f
    Sigma = Sigma_f - Kgain * C * Sigma_f;
    Zmat.col(tIndex) = mu;

    // 存储滤波后的估计结果
    mu_filtered[tIndex] = mu;
    Sigma_filtered[tIndex] = Sigma;
  }

  // 后向过程：Rauch-Tung-Striebel 平滑器
  mu_smooth = Zmat;
  Sigma_smooth = Sigma;
  for (tIndex = tTotal - 2; tIndex != (size_t)(-1); tIndex--) {
    // 计算平滑增益 J = Sigma_filtered[tIndex] * A^T * Sigma_f^{-1}
    mat J = Sigma_filtered[tIndex] * A.t() * inv(Sigma_f);

    // 更新平滑状态均值 mu_smooth
    mu_smooth.col(tIndex) = mu_filtered[tIndex] + J * (mu_smooth.col(tIndex + 1) - A * mu_filtered[tIndex]);
    // 更新平滑状态协方差 Sigma_smooth
    Sigma_smooth = Sigma_filtered[tIndex] + J * (Sigma_smooth - Sigma_f) * J.t();
    Zmat.col(tIndex) = mu_smooth.col(tIndex);
  }

  // 返回状态平滑估计矩阵 Zmat
  return Zmat;
}
```

### Rresult

```R
generateData = function(n, A, C, z0, sigmaQ, sigmaR) {
    # 定义函数 generateData，用于生成观测数据和实际的隐藏状态
    # 输入参数：
    # n - 观测的样本数
    # A - 状态转移矩阵
    # C - 观测矩阵
    # z0 - 初始状态向量
    # sigmaQ - 过程噪声标准差
    # sigmaR - 观测噪声标准差

    z = z0
    # 初始化实际隐藏状态矩阵 ZActual，大小为 4 x n
    ZActual = matrix(0, nrow = 4, ncol = n)
    # 初始化观测矩阵 Yobs，大小为 2 x n
    Yobs = matrix(0, nrow = 2, ncol = n) 
    
    # 模拟 n 个观测值
    for(i in 1:n) {
        # 根据状态转移方程更新隐藏状态 z = A * z + 过程噪声
        z = A %*% z + rnorm(4, sd = sigmaQ)
        # 计算观测值 Yobs = C * z + 观测噪声
        Yobs[, i] = C %*% z + rnorm(2, sd = sigmaR)
        # 存储当前时间步的隐藏状态
        ZActual[, i] = z
    }
    
    # 返回观测值和隐藏状态作为一个列表
    dataList = list(Yobs = Yobs, ZActual = ZActual)
    return(dataList)
}

drawPath = function(Yobs, ZActual, Zhat1, Zhat2) {
    # 定义函数 drawPath，用于绘制观测数据、实际隐藏状态和卡尔曼滤波、平滑的估计值
    # 输入参数：
    # Yobs - 观测矩阵
    # ZActual - 实际隐藏状态矩阵
    # Zhat1 - 卡尔曼滤波的状态估计
    # Zhat2 - 卡尔曼平滑的状态估计
    
    # 绘制观测数据和隐藏状态路径
    plot(t(Yobs), pch = 20, col = "grey", xlab = "", ylab = "")
    lines(t(ZActual[1:2, ]), col = "black", lwd = 2)
    lines(t(Zhat1[1:2, ]), col = "blue", lwd = 2)
    lines(t(Zhat2[1:2, ]), col = "red", lwd = 2)
}

# 参数设置
sigmaQ = 0.2  # 过程噪声标准差
sigmaR = 0.8  # 观测噪声标准差
n = 100       # 样本数量
delta = 0.5   # 步长
A = diag(1, nrow = 4, ncol = 4) # 状态转移矩阵
A[1, 3] = A[2, 4] = delta 
C = diag(1, nrow = 2, ncol = 4) # 观测矩阵
x0 = c(0, 0)  # 初始位置
v0 = c(1, 1)  # 初始速度
z0 = c(x0, v0) # 初始状态向量

# 去除注释符号以执行 kalmanSmoothing 函数
Rcpp::sourceCpp("homework2/homework2.3/kalmanFiltering.cpp")
Rcpp::sourceCpp("homework2/homework2.3/kalmanSmoothing.cpp")

# 生成观测数据
kdata = generateData(n, A, C, z0, sigmaQ, sigmaR)
Zhat1 = kalmanFiltering(kdata$Yobs, z0, A, C, sigmaQ, sigmaR)
Zhat2 = kalmanSmoothing(kdata$Yobs, z0, A, C, sigmaQ, sigmaR)

Yobs = kdata$Yobs
ZActual = kdata$ZActual

# 保存卡尔曼滤波和平滑结果图像
jpeg(filename = "~/homework2/homework2.3/Kalman_filtet_smoothing.png", width = 600, height = 600)

par(mar = c(2, 2, 2, 1), mfrow = c(2, 1))
plot(t(Yobs), pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Filter")
lines(t(ZActual[1:2, ]), col = "grey", lwd = 2)
lines(t(Zhat1[1:2, ]), col = "blue", lwd = 2)
plot(t(Yobs), pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Smoothing")
lines(t(ZActual[1:2, ]), col = "grey", lwd = 2)
lines(t(Zhat2[1:2, ]), col = "blue", lwd = 2)

dev.off()

# 显示卡尔曼滤波和平滑结果
par(mar = c(2, 2, 2, 1), mfrow = c(2, 1))
plot(t(Yobs), pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Filter")
lines(t(ZActual[1:2, ]), col = "grey", lwd = 2)
lines(t(Zhat1[1:2, ]), col = "blue", lwd = 2)
plot(t(Yobs), pch = 20, col = "grey", xlab = "", ylab = "", main = "Kalman Smoothing")
lines(t(ZActual[1:2, ]), col = "grey", lwd = 2)
lines(t(Zhat2[1:2, ]), col = "blue", lwd = 2)
```



## YauYau Filter

### 参数设置

```R
Dim <- 3
T <- 20
Dt <- 0.001
Dtau <- 5 * Dt
Nt <- as.integer(Dtau / Dt)
Ntau <- as.integer(T / Dtau)
NtNtau <- as.integer(T / Dt)

# Define functions f and h
f <- function(x) {
  return(c(cos(x[1]), cos(x[2]), cos(x[3]))) # R index starts from 1
}

h <- function(x) {
  return(c(x[1]^3, x[2]^3, x[3]^3)) # R index starts from 1
}
```

### 状态-观测模拟

```c++
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Simulate_State_Obser(double Dt, int Ntau, int NtNtau, Function f, Function h, int Dim) {
  mat x(NtNtau, Dim, fill::zeros);
  mat y_Dt(NtNtau, Dim, fill::zeros);
  mat y_tau(Ntau, Dim, fill::zeros);
  double sqrtdT = sqrt(Dt);
  
  // Set seed for reproducibility
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(1234);
  
  // Simulate x
  for (int t = 1; t < NtNtau; ++t) {
    rowvec x_prev = x.row(t - 1);
    NumericVector f_value = as<NumericVector>(f(x_prev));
    rowvec noise = randn<rowvec>(Dim);
    x.row(t) = x_prev + as<rowvec>(f_value) * Dt + sqrtdT * noise;
  }
  
  // Simulate y_Dt
  for (int t = 1; t < NtNtau; ++t) {
    rowvec x_prev = x.row(t);
    NumericVector h_value = as<NumericVector>(h(x_prev));
    rowvec noise = randn<rowvec>(Dim);
    y_Dt.row(t) = y_Dt.row(t - 1) + as<rowvec>(h_value) * Dt + sqrtdT * noise;
  }
  
  // Simulate y_tau
  for (int n = 1; n < Ntau; ++n) {
    int t = n * (NtNtau / Ntau);
    y_tau.row(n) = y_Dt.row(t);
  }
  
  return List::create(Named("x") = x,
                      Named("y_Dt") = y_Dt,
                      Named("y_tau") = y_tau);
}
```

```R
# 状态过程的图像
plot_State <- function(x) {
    # Convert matrices to data frames for easier plotting
    x_df <- as.data.frame(x)
    
    # Add time column
    x_df$time <- 1:nrow(x_df)
    
    # Reshape data frames for plotting
    x_melt <- reshape2::melt(x_df, id.vars = "time", variable.name = "Dimension", value.name = "Value")
    
    # Plot x for each dimension
    p1 <- ggplot(subset(x_melt, Dimension == "V1"), aes(x = time, y = Value)) +
      geom_line(color = "blue") +
      labs(title = "x_1", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    p2 <- ggplot(subset(x_melt, Dimension == "V2"), aes(x = time, y = Value)) +
      geom_line(color = "red") +
      labs(title = "x_2", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    p3 <- ggplot(subset(x_melt, Dimension == "V3"), aes(x = time, y = Value)) +
      geom_line(color = "green") +
      labs(title = "x_3", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Arrange the plots
    grid.arrange(p1, p2, p3, ncol = 3)
  }
```

```R
# 观测过程的图像
plot_Obser <- function(y_tau) {
    # Convert matrices to data frames for easier plotting
    y_tau_df <- as.data.frame(y_tau)
    
    # Add time column
    y_tau_df$time <- seq(1, nrow(y_tau_df))
    
    # Reshape data frames for plotting
    y_tau_melt <- reshape2::melt(y_tau_df, id.vars = "time", variable.name = "Dimension", value.name = "Value")
    
    # Plot y_tau for each dimension
    p4 <- ggplot(subset(y_tau_melt, Dimension == "V1"), aes(x = time, y = Value)) +
      geom_line(color = "blue") +
      labs(title = "y_1", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    p5 <- ggplot(subset(y_tau_melt, Dimension == "V2"), aes(x = time, y = Value)) +
      geom_line(color = "red") +
      labs(title = "y_2", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    p6 <- ggplot(subset(y_tau_melt, Dimension == "V3"), aes(x = time, y = Value)) +
      geom_line(color = "green") +
      labs(title = "y_3", x = "Time", y = "Value") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Arrange the plots
    grid.arrange(p4, p5, p6, ncol = 3)
  }
```

### PDE 离散

```C++
\\ 

\\ [[Rcpp::export]]
KolmogorovPDE_Discrete(){
\\

}
```

### PDE 求解

```C++
\\ 

\\ [[Rcpp::export]]
KolmogorovPDE_Solve(){
\\

}
```

### 状态估计

```C++
\\ 

\\ [[Rcpp::export]]
Estimation_State(){
\\
for (){
    \\
    \\
    for (){
        \\
        \\
    }
}

}
```

### 结果对比

```C++
\\ 结果对比

\\ [[Rcpp::export]]
Result_Compare(){
\\

}
```
