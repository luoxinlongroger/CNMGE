% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on July 1, 2022, in order to
% process the case without the analytical gradient of the objective
% function. 
% In this case, we use the automatic differentiation technique
% (dlgradient.m) to compute the analytical gradient of the objective
% function with the reverse mode in deep learning toolbox. 
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [opt_x,opt_f,CPU_time,Num_stp] = CNMGE(func,x0,ub,lb);
%
% This solver finds the global minimum of an unconstrained optimization 
% problems by using continuation Newton methods with deflation techniques 
% and quasi-genetic evolutionary algorithms. This program runs all programs 
% successively in desired way.
%
% This is the main subroutine of this software package. 
%
% Have a look at "Xin-long Luo, Hang Xiao, Sen Zhang
% Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% http//doi.org/10.21203/rs.3.rs-1102775/v1 or arXiv preprint available at
% http://arxiv.org/abs/2107.13864, and 
%"X.-L. Luo and H. Xiao, Generalized continuation Newton methods
% and the trust-region updating strategy for the underdetermined system,
% Journal of Scientific Computing, Vol. 88, 56 (2021), published online
% at http://doi.org/10.1007/s10915-021-01566-0, pp. 1-22, July 13, 2021.
% to see what this code does.
%
% The code can be available at
% "https://teacher.bupt.edu.cn/luoxinlong/zh_CN/zzcg/41406/list/index.htm"
% 
% Input:
%
% func: It denotes the objective function: R^{n}-->R. Since we use automatic
% differentiation to solve its gradient, it includes the objective function 
% and its gradient. If the output of func(x) only has one parameter, it only 
% compute the objective function value at x. Otherwise, it return the objective
% function value and its gradient at x. 
%
% x0: It denotes the initial point.
%
% lb: It denotes the lower bound of the variable (optional).
%
% ub: It denotes the upper bound of the variable (optional).
%
% Output:
%
% x_opt: It denotes the global optimal point of the objective function f(x), 
% which is found by CNMGE.
%
% f_opt: It denotes the global optimal value of the objective function f(x)
% at x_opt. 
%
% Num_stp: It denotes the number of found stationary points of f(x). 
%
% CPU_time: It denotes the computational time of CNMGE.
%
function [f_opt,x_opt,CPU_time,Num_stp] = CNMGE(func,x0,ub,lb)
tic;

if nargin < 3
    ub = inf;
    lb = -inf;
end

% Firstly, we call CNMDTM.m to find stationary points of the objective function 
% from multi-start points. The detailed implementations of this subroutine are refered to the
% descriptions of CNMDTM.m. 
tol = 1e-6; % the default tolerance. 

[sol_mat,success_flag] = CNMDTM(func,x0,ub,lb,tol);

if (success_flag == 0)
    disp('CNMDTM fails to find a stationary point of this problem.');   
   
    % Modify it on 2023/4/11 by Sen Zhang
    Num_stp = 0;
else
    [~,Num_stp] = size(sol_mat);
    fprintf('The number of stationary points found by CNMGE is %d\n',Num_stp);
end

% Then, we call QGE.m to approach the global optimal point closer.  
% The detailed implementations of this subroutine are refered to the
% descriptions of QGE.m. 

%modify on 2023/4/11
n = length(x0);
[f_opt, x_opt] = QGE(func, sol_mat, n, lb, ub);

% Get the computational time of CNMGE.
CPU_time = toc;

end

