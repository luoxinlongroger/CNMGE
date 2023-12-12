% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% The code can be available at
% "https://teacher.bupt.edu.cn/luoxinlong/zh_CN/zzcg/41406/list/index.htm"
%
% function J_new = jaco_deflation(J_pr,F,x,solx,n2,n1)
%
% This solver detetermines the global minmum of unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 
% and quasi-genetic ecolution. This program runs all programs successively 
% in desired way.
%
% This is the subroutine for calculating Jacobian matrix after using 
% deflation techneque of this software package. 
%
% Since J_k(x) with deflation techneque has the special structure, we can
% compute its sub-differential J_df of J_k(x) as follows:
%
%     J_df = cl/(||x-x_1||...||x-x_k||)*(J_pr + grad_func(x_k)*p(x)'),  (1)
%
% where the column vector function p(.) is defined by
%
%      p(x)' = -(sgn(x-x1)/||x-x_1|| + ... + sign(x-xk)/||x-x_k|| ),    (2)
%
% Have a look at title="Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% doi:10.21203/rs.3.rs-1102775/v1to, to see what this code does.
% 
% Input:
%
% J_pr: It denotes the the Jacobian matrix of the primitive problem at the
% current x_k point.
%
% grad_func: It denotes the gradient of the test problems.
%
% x_k: It denotes the current x_k.
%
% solx: It denotes the the set of local optimal solution obtained 
% by using CNMTD with the definition technique and switching updating
% technique.
%
% cl: It denotes the factor to prevent overflow.
%
% nl: It denotes the value of ||x-x_1||...||x-x_k||.
%
% Output: 
% J_df: It denotes the Hessian matrix of nonlinear equation with deflation
% technique.
%
function J_df = jaco_deflation(J_pr,func,x_k,solx,cl,nl)
% Get the dimensions of solution_mat where sm is the number of local
% optimal solutions obtained by CNMDT.
[~,Num_stp] = size(solx);
sum_dnorm = 0;

% F_k is the function value at x_k.
[~,Fx] = func(x_k);

% According to the special structure of J_df, we divide J_df into two
% parts:J_df1 and J_df2, we calculate them respectively.
% Calculate the value of J_df1.
J_df1 =J_pr*cl/nl;

% Calculate the value of J_df2.
for i = 1 : Num_stp 
    dnorm = sign(x_k-solx(:,i))/norm(x_k-solx(:,i),1);
    sum_dnorm = sum_dnorm + dnorm;
end
J_df2 = Fx*sum_dnorm'*cl*(-1/(nl));

% Merge J_df1 and J_df2 into J_df.
J_df = J_df1+J_df2;
end