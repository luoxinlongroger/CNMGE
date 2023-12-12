% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on August 25, 2022, in order to
% process the case without the given analytical gradient of the function. 
% In this case, we will compute the analytcial gradient of the objective
% function by the automatic differentiation technique (prob2struct.m) with the 
% reverse mode in the optimization toolbox. 
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [sol_x,successdt_flag,iter] = CNMDT(func,Num_init,sol_mat,x0,tol)
%
% This solver finds the global minimum of an unconstrained optimization 
% problem by using the continuation Newton methods with deflation techniques. 
%
% This subroutine is the continuation Newton methods[2] with the deflation 
% technique[1] and the adaptively updating Jacobian technique [3].
% The used techniques can be referenced as follows:
% 
% Assume that x_k is a zero point of G_{k-1}(x). Then, we construct another
% function G_k(.) via eliminating this zero point x_k of G_{k-1}(x) as
% fllows:
%
%                 G_k(x) = G_{k-1}(x)/(||x-x_k||)                    (1)
%
% where G_0(x) = g(x), g(x) is the gradient of f(x). Thus, from
% equation(1), we obtain
%
%               G_k(x) = g(x)/(||x-x_1||...||x-x_k||)                 (2)
%
% Have a look at title="Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% doi:10.21203/rs.3.rs-1102775/v1to, see what this code does.
%
% internal parameters:
%         itmax is the maximum number of iterations allowed.
%         epsilon is the convergence tolerance.
%         etag is the good approximation measurement of trust-region scheme.
%         etab is the bad approximation measurement of trust-region scheme. 
%         gammae is the enlarged factor of time-stepping size delta t.
%         gammar is the reduced factor of time-stepping size delta t.
%         rp0 is the regularization parameter.
%         Tol_Fk is the terminated tolerance of \|F_k - F_k1\|. 
% 
% Input:
%
% func: It denotes the objective function: R^{n}-->R. Since we use automatic
% differentiation to solve its gradient, it includes the objective function 
% and its gradient. If the output of func(x) only has one parameter, it only 
% compute the objective function value at x. Otherwise, it return the objective
% function value and its gradient at x. 
%
% sol_mat: It denotes the known stationary points of the objective function. 
%
% Num_init: It denotes the number of tried initial points. 
%
% x0: It denotes the initial point (optional).
%
% tol: It denotes the tolerance (optional). 
%
% Output: 
% 
% sol_x: It denotes the found stationary point of f(x). 
%
% successdt_flag: It denotes whether CNMDT finds a new statoinary
% point of f(x) or not. If successdt_flag be ture, it indicates that 
% a new stationary point is found by CNMDT. Otherwise, it indicates that 
% CNMDT fails to find a stationary point of f(x). 
%
% References:
% [1] Brown, K. M., Gearhart, W. B.: Deflation techniques for the 
% calculation of further solutions of a nonlinear system, Numer. Math. 16 
% (1971), 334-342.
%
% [2] Luo, X.-L., Xiao, H., Lv, J.-H.: Continuation Newton methods with the
% residual trust-region time stepping scheme for nonlinear equations, 
% Numer. Algorithms, Vol. 89, pp. 223-247, 2022, published online at 
% http://doi.org/10.1007/s11075-021-01112-x, 1-25, April 2, 2021.
%
% [3] Luo, X.-L., Xiao, H.: Generalized continuation Newton methods and the 
% trust-region updating strategy for the underdetermined system, Journal of 
% Scientific Computing, Vol. 88, article 56 (2021), published online at 
% https://doi.org/10.1007/s10915-021-01566-0, or available at 
% http://arxiv.org/abs/2103.05829, pp. 1-22, July 13, 2021.
%
% The code can be available at
% "https://teacher.bupt.edu.cn/luoxinlong/zh_CN/zzcg/41406/list/index.htm"
%
function [sol_x,successdt_flag,iter] = CNMDT(func,Num_init,sol_mat,x0,tol)
% Check the number of input parameters so that cnmdt can use different
% solution strategies.
%
if nargin < 5
    epsilon = 1.e-6;
else
    epsilon = tol;
end

% xm is the length of vector x_k.
[xm,Num_stp] = size(sol_mat);
if nargin < 4
    x0 = ones(xm,1); % Use the default initial point. 
end

if nargin < 3
    disp('It must input the objective function and known stationary points');
end

cm = xm; % cm is the length of gradient.

% Let initial point x_k equals to x0.
x_k = x0;

% Initialize the time step size delt_tk.
delt_tk = 1e-2;

% Initialize the delt_xk size.
delt_xk = 10;

% Initialize the number of iteration step.
iter = 0;

% the trust region parameters etag, etab.
etag = 0.25; % good approximation 
etab = 0.75; % bad approximation 

gammae = 2; % Enlarge delta t.
gammar = 0.5; % Reduce delta t.  

% Initialize the rho_k, rho_k is the ration of Actual reduction and
% predicted reduction.
rho_k = 0;

% Generate an identity matrix which dimension is cm * cm.
I = eye(cm,cm);

% Let regularization parameter rp0 equals to 1e-6.
rp0 = 1e-6;

% Initialize the cnt_Deltter.
cnt_Deltter = 0;

% The Termination tolerance on the function value.
Tol_Fk = 1e-6;


% Let the maximum number of iterations equal 200. 
itmax = 200;

% The following steps are the initialization step of the definition
% technology.
nl = 1; % nl means ||x-x_1||...||x-x_k||.
cl = 1; % cl is the factor to prevent overflow.
for i = 1 : Num_stp
    alph_k = norm(sol_mat(:,i),1);
    xdiffsol_n1 = norm(x_k-sol_mat(:,i),1);
    
    % If the norm of the ith solution is close to 0, multiplying it
    % directly into CL will invalidate cl. therefore, we set the norm
    % close to 0 to 1.
    if( alph_k < epsilon)
        alph_k = xm;
    end
    cl = cl*alph_k;
    nl = nl*xdiffsol_n1;
end

% When the inital point x0 equals the stationary point, we cannot use the 
% deflation technique and switch to another initial poit, then we call CNMDT 
% again, 2022 Oct. 17, Xin-long Luo revised it. 
if nl < 1.e-10
    sol_x = x_k;
    successdt_flag = 0;
    return;
end


% Res_k is set to the residue with deflation technique.
[~,gxk] = func(x_k);
F_k = gxk*cl/nl;
Res_k = norm(F_k,inf);

% Use a while loop to solve nonlinear equation cnsteq.
while(Res_k>epsilon)      
    % Due to the reason of deflation technique, the behavior of the
    % problem will gradually become worse. Therefore, the judgment
    % conditions of the adaptively updating Jacobian matrix need to be more
    % strict.
    if( (abs(1-rho_k) > etag)|| delt_xk > 1 || (Num_init > 2 && Num_stp < 5))
        % we call jacobianfun.m to calculate the Jacobian matrix of the
        % primitive nonlinear equation.
        J_pr = jacobianfun(func,x_k);
        % we call jaco_deflation.m to calculate the Jacobian matrix of 
        % nonlinear equation with deflation technique.
        J_xk = jaco_deflation(J_pr,func,x_k,sol_mat,cl,nl);
    end
    
    % In order to increase the computational stability, we add the
    % regularization alph0*I to J_xk. We use the explicit continuation 
    % Newton methods to solve s_k.
    s_kp = -(J_xk+rp0*I)\F_k;%
    s_k = (delt_tk/(1+delt_tk))*s_kp;
    
    % Update the nl.
    nl = 1;
    for i = 1 : Num_stp
        xdiffsol_n1 = norm(x_k+s_k-sol_mat(:,i),1);
        nl = nl*xdiffsol_n1;
    end
    
    % The following steps are the trust-region scheme for adjusting
    % time-stepping size 
    x_kp = x_k + s_k;
    [~,gxkp] = func(x_kp);
    F_kp = gxkp*cl/nl;

    normF_k = norm(F_k,2);
    normF_kp = norm(F_kp,2);
    % Compute the predicted reduction of the residue.
    Pred_k = normF_k - normF_kp; 

    % Compute the actual reduction of the residue.
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
    % If the gradient of the problem remains unchanged, the calculation
    % is stopped.
    delt_Fk = norm(F_kp-F_k);
    if(delt_Fk<Tol_Fk)
        cnt_Deltter = cnt_Deltter + 1;
        if(cnt_Deltter> 10)
            break;
        end
    else
        cnt_Deltter = 0;
    end
    
    % Update the x_k, residual and the delt_xk of the problem.
    delt_xk = norm(s_k,inf);
    x_k = x_kp;
    F_k = F_kp;   
    Res_k = norm(F_k,inf);   
    
    % If the time step size is too small, stop CNMDT.
    if(delt_tk < 1e-7)
        break;
    end
    
    % If iter exceeds the upper limit, stop CNMDT.
    iter = iter + 1;
    if(iter >= itmax)
        break;
    end
end
sol_x = x_k;
if (Res_k < epsilon)
    successdt_flag = 1;
else
    successdt_flag = 0;
end
end

