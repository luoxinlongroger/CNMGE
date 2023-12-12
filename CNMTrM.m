% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on August 25, 2022, in order to
% process the case without the given analytical gradient of the function. 
% In this case, we will compute the analytcial gradient of the objective
% function by the automatic differentiation technique (prob2struct.m) with the 
% reverse mode in the optimization toolbox. 
% This code is written by Xin-long Luo on Sep. 8, 2022. 
% CNMTr finds a stationary point of f(x). 
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [x_opt,f_opt,success_flag,iter] = CNMTrM(f,x0,tol,ub,lb)
%
% This solver finds the global minimum of an unconstrained optimization 
% problem by using the continuation Newton methods from multi-start points.
%
% This subroutine is the continuation Newton methods[1],[2] with the 
% trust-region updating strategy to find multiple stationary points of f(x)
% from the default multi-start points. 
%
% Have a look at title="Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% doi:10.21203/rs.3.rs-1102775/v1to, see what this code does.
%
% internal parameters:
%         maxit is the maximum number of iterations allowed.
%         epsilon is the convergence tolerance.
%         etag is the good approximation measurement of the trust-region scheme.
%         etab is the bad approximation measurement of the trust-region scheme. 
%         gammae is the enlarged factor of time step delta_t.
%         gammar is the reduced factor of time step  delta_t.
%         rp0 is the regularization parameter.
%         Tol_Fk is the terminated tolerance of \|F_k - F_k1\|, where F_k is the 
%         gradient of f(x_k), and F_k1 is the gradient of f(x_k1).
% 
% Input:
%
% f: It denotes the objective function: R^{n}-->R. Since we use automatic
% differentiation to solve its gradient, it includes the objective function 
% and its gradient. If the output of f(x) only has one parameter, it only 
% compute the objective function value at x. Otherwise, it return the objective
% function value and its gradient at x. 
%
% x0: It denotes the given initial point.
% 
% tol: It denotes the tolerance (optional). 
%
% J: J(x) denotes the Jacobian matrix of \nabla f(x), i.e. 
% the Hessian matrix of the objective function f(x) (optional).
% If it does not given, we use the finite difference to compute it. 
%
% Output: 
%
% x_opt: It denotes a optimal approximation found by CNMTrM. 
%
% f_opt: It denotes the objective function value f(x_opt) at 
% the optimal approximation x_opt found by the CNMTrM subroutine. 
%
% success_flag: It indicates whether a stationary point is found by CNMTr
% or not. If success_flag is 1, it indicates a stationary point is found by
% CNMTr successfully. Otherwise, it indicates that CNMTr fails to find a
% stationary point.
%
%
% References:
% [1] Luo, X.-L., Xiao, H., Lv, J.-H.: Continuation Newton methods with the
% residual trust-region time stepping scheme for nonlinear equations, 
% Numer. Algorithms, Vol. 89, pp. 223-247, 2022, published online at 
% http://doi.org/10.1007/s11075-021-01112-x, 1-25, April 2, 2021.
%
% [2] Luo, X.-L., Xiao, H.: Generalized continuation Newton methods and the 
% trust-region updating strategy for the underdetermined system, Journal of 
% Scientific Computing, Vol. 88, article 56 (2021), published online at 
% https://doi.org/10.1007/s10915-021-01566-0, or available at 
% http://arxiv.org/abs/2103.05829, pp. 1-22, July 13, 2021.
%
% [3] X.-L. Luo and H. Xiao: Continuation Newton methods with deflation
% techniques and quasi-genetic evolution for global optimization problems,
% arXiv preprint available at http://arxiv.org/abs/2107.13864, 
% or Research Square preprint available
% at https://doi.org/10.21203/rs.3.rs-1102775/v1,  July 30, 2021.
%
% The code can be available at
% "https://teacher.bupt.edu.cn/luoxinlong/zh_CN/zzcg/41406/list/index.htm"
%
function [x_opt,f_opt,successmp_flag] = CNMTrM(f,x0,tol,ub,lb)
% Let the default initial point x_k equal x0.
if nargin < 3
    tol = 1.e-6;
end
%
if nargin < 2
    disp('It must input the objective function and its initial point');
end


successmp_flag = 0; % It denots that the continuation Newton (CNMTrM) fails to 
                    % find a stationary point of f(x) from multi-start points. 

N = 6; % the number of initial points.

n = length(x0);

solg_mat = zeros(n,N); % It keeps all found stationary points of f(x). 
fopt_vec = zeros(N,1); % It denotes the objective function value at the stationary 
                   % point x_opt. 
Num_stp = 0; % It denotes the number of the staionary points. 

for Num_init = 1 : N
    x0_p = Switch_Initial_Point(x0, Num_init);
    [x_spt,f_spt,success_flag] = CNMTr(f,x0_p,tol);    
 
    % If the continuation Newton method (CNMTr) finds a stationary point of
    % f(x) successfully, we let Num_init increase one and keep this stationary
    % point. 
    if (success_flag == 1)
        Num_stp = Num_stp + 1; 
        solg_mat(:,Num_stp) = x_spt;
        fopt_vec(Num_stp) = f_spt; 
    end
end

i = 1;
while(i <= N)
    % If the maximum value of the solution is greater than the upper 
    % boundary or the minimum value is lower than the lower boundary, 
    % the solution is removed from sol_mat_seeds.
    if(i<=Num_stp)
        % Get the maximum value in the ith solution.
        max_x = max(solg_mat(:,i));

        % Get the miximum value in the ith solution.
        min_x = min(solg_mat(:,i));
        if(max_x>ub || min_x<lb)
            solg_mat(:,i) = [];
            Num_stp = Num_stp-1;
            fopt_vec(i,:) = [];
            i = i-1;
            N = N-1;
        end
    else 
        solg_mat(:,i) = [];
        fopt_vec(i,:) = [];
        i=i-1;
        N=N-1;
    end
    
    % If the maximum and minimum values of the solution satisfy the upper 
    % and lower boundaries, the solution is retained.
    i = i + 1;
end

Num_stp = length(fopt_vec);
if (Num_stp < 1)
   x_opt = x0;
   f_opt = f(x0);    
else
    successmp_flag = 1; % It denots that the continuation Newton (CNMTrM)  
                        % finds at least a stationary point of f(x) from 
                        % multi-start points. 
%     f_opt = fopt_vec(1);
%     x_opt = solg_mat(:,1);
%     for i = 2:Num_stp
%         if (f_opt > fopt_vec(i))
%             f_opt = fopt_vec(i);
%             x_opt = solg_mat(:,i);
%         end
%     end  

    [f_opt,idx] = min(fopt_vec);
    x_opt = solg_mat(:,idx);
end 
end

