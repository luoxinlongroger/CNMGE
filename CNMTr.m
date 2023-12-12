% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on August 25, 2022, in order to
% process the case without the given analytical gradient of the function. 
% In this case, we will compute the analytcial gradient of the objective
% function by the automatic differentiation technique (prob2struct.m) with the 
% reverse mode in the optimization toolbox. 
% This code is written by Xin-long Luo on Sep. 8, 2022, according to the 
% revised paper. Namely, he divides CNMDT into three functions: CNMTR, CNMDT 
% and CNMDTM. 
% CNMTR finds a stationary point of f(x). 
% CNMDT finds multiple stationary points of f(x) from an initial point. 
% CNMDTM finds multiple stationary points of f(x) from multiple initial
% points.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [x_spt,f_spt,success_flag,iter] = CNMTr(f,x0,tol)
%
% This solver finds the global minimum of an unconstrained optimization 
% problem by using the continuation Newton methods. 
%
% This subroutine is the continuation Newton methods[1],[2] with the 
% trust-region updating strategy to find a stationary point of f(x). 
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
% x_spt: It denotes a stationary point found by CNMTr. 
%
% f_spt: It denotes the objective function value f(x_spt) at 
% the stationary point x_spt found by the CNMTr subroutine. 
%
% success_flag: It indicates whether a stationary point is found by CNMTr
% or not. If success_flag is 1, it indicates a stationary point is found by
% CNMTr successfully. Otherwise, it indicates that CNMTr fails to find a
% stationary point.
%
% iter: It denotes iterations of CNMTr.
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
function [x_spt,f_spt,success_flag,iter] = CNMTr(f,x0,tol)
%
if nargin < 3
    tol = 1.e-6;
end
%
if nargin < 2
    disp('It must input the objective function and give an initial point');
end

epsilon = tol;

% Let the default initial point x_k equal x0.
if nargin < 2
    disp('It must input the objective function and its initial point');
else
    x_k = x0;
end

% xm is the length of vector x_k.
[xm,~] = size(x_k);

cm = xm; % cm is the length of gradient.

% Initialize the time step delt_tk.
delt_tk = 1e-2;

% Initialize the number of iteration step.
iter = 0;

% The default value of success_flag is 0 and it indicates that CNMTr 
% fails to find a stationary point of f(x). 
success_flag = 0; 

% trust region parameter etag, etab.
etag = 0.25; % good approximation 
etab = 0.75; % bad approximation 

% the enlarging factor of time step gamma1, gamma2, gamma3
gammae = 2; % enlarge delta t 
gammar = 0.5; % reduce delta t 

% Initialize the rho_k. rho_k is the ratio of the actual reduction to 
% the predicted reduction.
rho_k = 0;

% Generate an identity matrix which dimension is cm * cm.
I = eye(cm,cm);

% Set the default regularization parameter rp0 and let it equal 1e-6.
rp0 = 1e-6;

% Initialize the cnt_Deltter.
cnt_Deltter = 0;

% The Termination tolerance of the function value.
Tol_Fk = 1e-6;

% Let the default maximum number of iterations equal 200.
maxit = 200;

[~,F_k] = f(x_k);
Res_k = norm(F_k,inf); % Compute the infinite norm of the gradient F_k.

% Use a while loop to find a stationary point of f(x). 
while(Res_k > epsilon) 
    % If (abs(1-rho_k) > etag), which means F_k+J_k*s_k approximates
    % F(x_k + s_k) well. We let J_{k+1} = J_k. Otherwise, we let J_{k+1} 
    % J(x_{k+1}). 
    if( (abs(1-rho_k) > etag))
        % we call jacobianfun.m to calculate the Jacobian matrix of the 
        % gradient g(x), i.e. the Hessian matrix of the objective function
        % f(x). 
        J_xk = jacobianfun(f,x_k);
    end
            
    % In order to increase the computational stability, we add the
    % regularization alph0*I to J_xk. We use the explicit continuation 
    % Newton methods to solve S_k.        
    s_kp = -(J_xk + rp0*I)\F_k;
    s_k = (delt_tk/(1+delt_tk))*s_kp;
    
    % The following steps are the trust-region scheme for adjusting
    % time step adaptively.  
    x_kp = x_k + s_k;
    % Compute the gradient of the objective at the predicted point
    % x_kp. 
    [~,F_kp] = f(x_kp);
    
    normF_k = norm(F_k,2);
    normF_kp = norm(F_kp,2);
    % Compute the predicted reduction of the gradient length.
    Pred_k = normF_k - normF_kp; 

    % Compute the actual reduction of the gradient length.
    Ared_k = (delt_tk/(1+delt_tk))*normF_k;

    % Compute the ratio of the actual reduction to the 
    % predicted reduction.
    rho_k = abs(Pred_k/Ared_k);  

    % Adjust the time step by the trust-region updating strategy 
    % adaptively. 
    if(abs(1-rho_k) > etab)
        % It is a bad approximation and reduce the time step, 
        % and retains the iteration step x_k and its gradient F_k.
        delt_tk = gammar*delt_tk;       % Reduce the time step. 
    end
    if(abs(1-rho_k) <= etag)
        % It is a good approximation. Enlarge the time step  
        % and accept the trial step.
        delt_tk = gammae*delt_tk;       % Enlarge the time step. 
    end
    
    % If the gradient of the objective function remains unchanged, it stops.     
    delt_Fk = norm(F_kp - F_k);
    if(delt_Fk < Tol_Fk)
        cnt_Deltter = cnt_Deltter + 1;
        if(cnt_Deltter > 10)
            break;
        end
    else
        cnt_Deltter = 0;
    end   
    
    x_k = x_kp; % Update x_k. 
    F_k = F_kp;
    Res_k = norm(F_k,inf); % Compute the infinite norm of the gradient F_k.

    % If iter exceeds the upper limit, stop CNMTr.
    if (iter > maxit)
        break;
    else
        iter = iter + 1;
    end        
end
if (Res_k <= epsilon)
    success_flag = 1; 
end
x_spt = x_k;
f_spt = f(x_spt);
end

