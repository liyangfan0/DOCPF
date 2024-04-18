function [mu, Sigma]=kalman(mu,Sigma,obs,Q,R,A,C)
%mu均值，Sigma协方差矩阵，obs观测值，Q过程噪声，R观测噪声
%A是由运动模型确定的，C是由观测模型确定的
[mu_bar, Sigma_bar] = kalmanPredict(mu, Sigma, A, Q);
[mu, Sigma] = kalmanUpdate(mu_bar, Sigma_bar, obs, C, R);
