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
% function [solg_mat,Num_stp,success_flag,x_opt,f_opt] = CNMDTM(func,x0,ub,lb,tol)
%
% This solver finds the global minimum of an unconstrained optimization 
% problem by using the continuation Newton methods with deflation techniques 
% from multi-start points. 
%
% This subroutine is the continuation Newton methods[2] with the deflation 
% technique [1] and the adaptively updating Jacobian technique[3].
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
% Input:
%
% func: It denotes the objective function: R^{n}-->R. Since we use automatic
% differentiation to solve its gradient, it includes the objective function 
% and its gradient. If the output of func(x) only has one parameter, it only 
% compute the objective function value at x. Otherwise, it return the objective
% function value and its gradient at x. 
%
% x0: It denotes the initial point (optional).
%
% tol: It denotes the tolerance (optional). 
%
% lb: It denotes the lower bound of the variable.
%
% ub: It denotes the upper bound of the variable.
%
% Output: 
% 
% sol_mat: It denotes all found stationary points of the objective function. 
%
% success_flag: It denotes whether CNMDTM finds at least a new statoinary
% point of f(x) or not. If success_flag be ture, it indicates that at
% least a new stationary point is found by CNMDTM. Otherwise, it indicates that 
% CNMDTM fails to find a stationary point of f(x). 
%
% x_opt: It denotes an optimal approximation point computed by CNMDTM. 
%
% f_opt: It denotes the objective function value at the optimal 
% approximation x_opt computed by the CNMDTM subroutine. 
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
function [sol_mat,success_flag,x_opt,f_opt] = CNMDTM(func,x0,ub,lb,tol)
%
if nargin < 5
    tol = 1.e-6;
end
%
if nargin < 4
    disp('It must input the objective function and give an initial point');
end

% It gives an default value, which means that the continuation Newton fails 
% to find a stationary point of f(x) from multi-start points. 
success_flag = 0; 

Num_stp = 0;  % Initialize the number of stationary points. 

Num_init = 6; % the number of initial points. 

% Firstly, we call CNMTr.m to find a stationary point of the objective function.
% The detailed implementation of this subroutine is refered to the descriptions
% of CNMTr.m. 

for No_init = 1 : Num_init
    x0_p = Switch_Initial_Point(x0, No_init);
    [x_spt,f_spt,success_flag,~] = CNMTr(func,x0_p,tol);    
 
    % If the continuation Newton method (CNMTr) finds a stationary point of
    % f(x) successfully, we let No_init increase one and keep this stationary
    % point. 
    if (success_flag == 1)
        Num_stp = Num_stp + 1; 
        sol_mat = x_spt;
        fopt_vec = f_spt;         
        break;    
    end
end

% We use CNMDT and a while loop to implement the deflation technique to find 
% stationary points of the objective function as many as possible from 
% multi-start points. After the maximum number of initial points are tried, 
% the loop is terminated.
if (success_flag == 1)    
    No_init = 1;    
    while(No_init <= Num_init)
        % Call CNMTD to find a new stationary point of the objective
        % function, which is different from the found stationary points.   
        x0_p = Switch_Initial_Point(x0, No_init);
        [sol_x,successdt_flag,~] = CNMDT(func,No_init,sol_mat,x0_p);
         
        % If CNMDT finds a new stationary point and it does not equal x0_p, 
        % we store this found stationary point in sol_mat and continue calling 
        % CNMDT to find the new stationary point from the same initial point. 
        % Otherwise, we switch a new initial point and call CNMDT to find 
        % a new stationary point.  
        if(successdt_flag == 1)             
            % we store the stationary point into sol_mat.
            if(norm(sol_x-x0_p,1) < tol)
                % Switch to the new initial point.  
                No_init = No_init + 1;
            end
            sol_mat_new = [sol_mat,sol_x];             
            sol_mat = sol_mat_new;  
            func_x = func(sol_x);
            fopt_vec_new = [fopt_vec,func_x];
            fopt_vec = fopt_vec_new;
            % The number of found stationary points increases by one.
            Num_stp = Num_stp +1; 
        else
            No_init = No_init + 1;
        end        
    end

%modify on 2023/4/11
else
    sol_mat = [];
%     return;
end

i=1;
while(i<=Num_stp)
    % Get the maximum value in the ith solution.
    max_x = max(sol_mat(:,i));

    % Get the miximum value in the ith solution.
    min_x = min(sol_mat(:,i));
    if(max_x > ub || min_x < lb)
        sol_mat(:,i) = [];
        Num_stp = Num_stp-1;
        fopt_vec(:,i) = [];
        i = i-1;
    end

    i = i+1;
end


if nargout > 3
    if (Num_stp < 1)
       x_opt = x0;
       f_opt = func(x0);   
       success_flag = 0;
    else
        success_flag = 1; % It means that the continuation Newton (CNMDTM)
        % finds at least a stationary point of f(x) from multi-start points.
        fopt_out_vec = fopt_vec(1:Num_stp);
        [f_opt,index_min] = min(fopt_out_vec);
         x_opt = sol_mat(:,index_min); 
    end
end
end



