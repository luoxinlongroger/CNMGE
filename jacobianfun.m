% This code is developed by Xinlong Luo. September, 2020.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function Mat_jacobian = jacobianfun(F,x)
%
% This solver detetermines the global minmum of unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 
% and quasi-genetic ecolution. This program runs all programs successively 
% in desired way.
%
% This is the subroutine for calculating Jacobian matrix of this software
% package. Use the vector operator of MATLAB to improve the computational
% efficiency.
%
% Have a look at "Xin-long Luo, Hang Xiao, Sen Zhang
% Continuation Newton methods with deflation techniques and quasi-genetic 
% evolution for global optimization problems", 
% https://doi.org/10.21203/rs.3.rs-1102775/v1, or 
% arXiv preprint available at http://arxiv.org/abs/2107.13864, 
% to see what this code does.
%
% The code can be available at
% "https://teacher.bupt.edu.cn/luoxinlong/zh_CN/zzcg/41406/list/index.htm"
% 
% Input:
%
% grad_func: It denotes the gradient of the test problems.
%
% x_k: It denotes the current x_k.
%
% Output: 
% Mat_jacobian: It denotes the Hessian matrix of nonlinear equation 
%
function Mat_jacobian = jacobianfun(func,x_k)
delta_x = 2.0e-8; 

% n is the length of vector x. 
n = length(x_k);

% Fx is the function value of grad_func at x_k.
[~,Fx] = func(x_k);

% m is the dimension of function Fx.
m = length(Fx);  

% Initialize the Jacobian matrix.
Mat_jacobian = zeros(m,n);

% Use for loop to calculating the Jacobian matrix.
for i = 1:n    
    incr_x = x_k;
    incr_x(i) = incr_x(i) + delta_x;
    [~,Fincr_x] = func(incr_x);
    Mat_jacobian(:,i) = (Fincr_x - Fx)/delta_x;
end
end
