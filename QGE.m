% This code is developed by Hang Xiao and Xinlong Luo. October, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [f_opt, x_opt] = QGE(func,sol_mat,lb,ub)
%
% This subroutine finds the global minimum of an unconstrained optimization 
% problem by quasi-genetic evolutionary algorithms and continuation Newton
% methods. 
%
% We use the idea of memetic algorithms to approach the global minimum further. 
% Namely, we use the heuristic crossover of the genetic algorithm to evolve
% for approaching the global minimum of the objective function from its local 
% optimal solutions, which are solved by CNMDTM.
%
% Have a look at "Xin-long Luo, Hang Xiao, 
% Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% http//doi.org/10.21203/rs.3.rs-1102775/v1 or arXiv preprint available at
% http://arxiv.org/abs/2107.13864, to see what this code does.
% 
% Input:
%
% sol_mat: It denotes the set of found stationary points by the subroutine 
% CNMTDM, which is the continuation Newton method with the deflation technique 
% from multi-start points. 
%
% func: It denotes the objective function. 
%
% lb: It denotes the lower bound of the variable.
%
% ub: It denotes the upper bound of the variable.
%
% Output: 
%
% f_opt: It denotes the global approximation minimum of the 
% objective function.
%
% x_opt: It denotes the global optimal approximation solution of the
% optimization problem.
%
function [f_opt, x_opt] = QGE(func, sol_mat, n, lb, ub)

% Obtain the dimensions of sol_mat. n is the dimensions of the problem and 
% Num_stp is the number of found statoinary points by CNMDTM.

% [n,Num_stp] = size(sol_mat);

%modify on 2023/4/11
[~,Num_stp] = size(sol_mat);

% If the number of known stationary points of f(x) is less than L, we add the 
% SUB_MAT to the solg_mat. There are L determistic seeds such as
% {zeros(n,1), 0.1*SUB_MAT, SUB_MAT, 10*SUN_MAT, 100*SUB_MAT,
% 1000*SUB_MAT}, where the SUB_MAT = {ones(n,1), -ones(n,1),
% [ones(n/2,1);-ones(n/2,1)], [-ones(n/2,1);ones(n/2,1)]}.

sub_seed0 = zeros(n,1);
sub_seed1 = ones(n,1);
sub_seed2 = -ones(n,1);
sub_seed3 = sub_seed1;
sub_seed3(n/2+1:end) = -sub_seed1(n/2+1:end);
sub_seed4 = sub_seed2;
sub_seed4(n/2+1:end) = -sub_seed2(n/2+1:end);
SUB_MAT = [sub_seed1 sub_seed2 sub_seed3 sub_seed4];
SUB_MAT = [sub_seed0 0.1*SUB_MAT SUB_MAT 10*SUB_MAT 100*SUB_MAT 1000*SUB_MAT];
[~,L] = size(SUB_MAT); % Add L initial candidate seeds. 
% sol_mat_seeds = [sol_mat SUB_MAT];
% Num_seeds = Num_stp + L; 

%modify on 2023/4/11
if(Num_stp ~= 0)
    sol_mat_seeds = [sol_mat SUB_MAT];
else 
    sol_mat_seeds = SUB_MAT;
end
Num_seeds = Num_stp + L; 

% Use while loop to remove the solutions that do not satisfy the upper and 
% lower bounds of x in the solution_mat.
% Initialize the counter i.
i = 1;
while(i <= Num_seeds)
    % Get the maximum value in the ith solution.
    max_x = max(sol_mat_seeds(:,i));
    
    % Get the miximum value in the ith solution.
    min_x = min(sol_mat_seeds(:,i));
    
    % If the maximum value of the solution is greater than the upper 
    % boundary or the minimum value is lower than the lower boundary, 
    % the solution is removed from sol_mat_seeds.
    if(max_x>ub || min_x<lb)
        sol_mat_seeds(:,i) = [];
        Num_seeds = Num_seeds-1;
        i = i-1;
    end
    
    % If the maximum and minimum values of the solution satisfy the upper 
    % and lower boundaries, the solution is retained.
    i = i + 1;
end
% If the population number is greater than the number of local optimal 
% solutions, let the population number equal to the number of local optimal 
% solutions.
if (L > Num_seeds)
    L = Num_seeds;
end

% Initializes the number of offspring populations.
solution_num = L*(L-1)/2;

% Initializes the f_solution_vec used to record the function value.
f_solution_vec = zeros(Num_seeds,1);

% Store the function values of the local optimal solutions in 
% f_solution_vec.
for i = 1 : Num_seeds
    f_solution_vec(i) = func(sol_mat_seeds(:,i));
end

% In order to make the effect of quasi genetic algorithms more outstanding, 
% we choose the first L minimizers of function as the initial population.
[f_solution_vec_new,idx] = sort(f_solution_vec);

% Get the first L minimizers of function values. 
f_solution_min_vec = f_solution_vec_new(1:L);

% Rearrange the sequence of the found stationary points according to the 
% ascendent index of their function values.
sol_mat_seeds_new = sol_mat_seeds(:,idx);
% sol_mat_seeds_new = sol_mat_seeds(:,);

% Select the first L minimum points. 
solution_min_mat = sol_mat_seeds_new(:,1:L);

% Initializes the number of iterations.
loop_iter = 1;

Num_Evo = 20; % the evolutionary generations. 

% Initializes the evolutionary solution mat.
e_solution_mat = zeros(n,solution_num);
    
% Initializes the evolutionary function value vector.
f_value_vec = zeros(solution_num,1);

% Use while loop to realize the cross evolution process of process of 
% quasi genetic algorithm and try to find the global optimal solution.
while(loop_iter < Num_Evo)   
    
    % Initializes the index of offspring.
    idx_ofs = 1;
    
    % Use ywo for loops to realize cross evolution between seeds and 
    % generate offspring.
    for i = 1 : L
        for j = i+1 : L
            % Get a new offspring by adding two different seeds and 
            % dividing by 2.
            new_offspring = ...
            (solution_min_mat(:,i)+solution_min_mat(:,j))/2;

            % Save the new offspring in the e_solution_mat.
            e_solution_mat(:,idx_ofs) = new_offspring;
            
            % Save the new offspring value in the f_value_vec.
            f_value_vec(idx_ofs) = func(new_offspring);
            
            % The number of idx_ofs increases by one.
            idx_ofs = idx_ofs + 1;
        end
    end
    % Merge the parents and the generated offspring into solution_min_mat.
    solution_min_mat_new = [solution_min_mat,e_solution_mat];
    solution_min_mat = solution_min_mat_new;
    
    % Merge the function values of parents and the function values of 
    % offspring into f_solution_min_vec.
    f_solution_min_vec_new = [f_solution_min_vec;f_value_vec];
    f_solution_min_vec = f_solution_min_vec_new;
    
    % We combine the local optimal solutions with the corresponding 
    % function values to reduce the amount of subsequent operations.
    % total_solution = [solution_min_mat;f_solution_min_vec];
    
    % Use the sort function to get the index of the ascendent function values.    
    [f_solution_min_vec_new,idx] = sort(f_solution_min_vec);

     % Get the first L minimizers of function. 
    f_solution_min_vec = f_solution_min_vec_new(1:L);
    
    % Rearrange the sequence of the local optimal solutions according to 
    % the index obtained in the previous step.
    solution_min_mat_new = solution_min_mat(:,idx); 

    % Get the first L minimum points.     
    solution_min_mat = solution_min_mat_new(:,1:L); 

    % The number of loop_iter increases by one.
    loop_iter = loop_iter + 1;
end
% Initializes the evolutionary result matrix.
solution_result = zeros(n,L);

% Initializes the evolutionary function value vector.
func_value_result = zeros(L,1);

% We refine the evolutionary results by using CNMTr from them. 
tol = 1.e-6;
for i = 1 : L    
    x0 = solution_min_mat(:,i);
    [x_spt,f_spt,~,~] = CNMTr(func,x0,tol);   
    solution_result(:,i) = x_spt;
    func_value_result(i) = f_spt;
end   

% Find the index of the minimum value in the func_value_result.
[f_opt,minf_index] = min(func_value_result);

% Set the smallest function value in the func_value_result as the 
% global_function_value.
x_opt = solution_result(:,minf_index);

% If the new computed global minimum is greater than that of genetic evolution, 
% the global approximation minimum is set by the minimum solution of genetic 
% evolution.
if(f_opt > f_solution_min_vec(1))
    f_opt = f_solution_min_vec(1);
    x_opt = solution_min_mat(:,1);
end

end

