% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on Augus 14-22, 2022, in order to
% process the case with the analytical gradient of the objective
% function.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% This solver finds the global minimum of an unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 
% and quasi-genetic evolutionary algorithms. This program runs all programs 
% successively in desired way.
%
% This is an entrance subroutine when we use this software to find the
% global optimal solution of an unconstrained optimization problem.
%
% Have a look at "Xin-long Luo, Hang Xiao, Sen Zhang, 
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
% There is no input parameter. However, any of the 68 test problems 
% provided can be selected for testing, and each test problem has 6 
% parameters:
% n: it denotes the dimension of the test case.
% 
% x0: it denotes the initial point of CNMGE.
%
% func: It denotes the objective function: R^{n}-->R. It includes the 
% objective function and its analytical gradient. If the output of 
% func(x) only has one parameter, it only compute the objective function 
% value at x. Otherwise, it return the objective function value and 
% its gradient at x. 
%
% lb: it denotes the lower boundary of the variable.
%
% ub: it denotes the upper boundary of the variable.
%
% Output: 
% There is no output parameter. We output the global optimal solution 
% of the problem, the number of local optimal solutions found by CNMGE
% and the computational time.
%
% We call gcp('nocreate') to get current parallel pool.
% poolobj = gcp('nocreate');% If no pool,  create new one.
% if isempty(poolobj)
%     poolsize = 0;
%     CoreNum=4; 
%     parpool(CoreNum);
% else
%     poolsize = poolobj.NumWorkers;
%     disp('Already initialized');
% end


clear;

test = input('\n Please input the number of test problem:');
switch test
    case 1 
        % molecular energy problem 
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function 
        func = @molecular_energy_problem;

        lb = -inf;
        ub = inf;
        
    case 2  
        % ackley_function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function 
        func = @ackley_function;

        lb = -inf;
        ub = inf;
        
    case 3
        % levy_function
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function 
        func = @levy_function;

        lb = -inf;
        ub = inf;
        
    case 4  
        % schwefel_function
        n = 1000;
        x0 = 200*ones(n,1);

        % Define the objective function 
        func = @schwefel_function;

        lb = -500;
        ub = 500;
        
    case 5
        % rastrigin_function 
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @rastrigin_function;

        lb = -inf;
        ub = inf;
        
    case 6
        % styblinski_tang_function 
        n = 1000;
        x0 = ones(n,1);

         % Define the objective function
        func = @styblinski_tang_function;

        lb = -inf;
        ub = inf;
        
    case 7
         % trid function
         n = 1000;
         x0 = ones(n,1);

         % Define the objective function
         func = @trid_function;

         lb = -inf;
         ub = inf;
         
    case 8
        % sum_squares_function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @sum_squares_function;

        lb = -inf;
        ub = inf;
        
     case 9
        % sphere_function 
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @sphere_function;

        lb = -inf;
        ub = inf;
        
     case 10
        % rotated_hyper_ellipsoid_function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @rotated_hyper_ellipsoid_function;

        lb = -inf;
        ub = inf;
        
     case 11
        % zakharov_function 
        n = 1000;
        x0 = (5e-5)*ones(n,1);

        % Define the objective function
        func = @zakharov_function;

        lb = -inf;
        ub = inf;
        
     case 12
        % dixon_price_function 
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @dixon_price_function;

        lb = -inf;
        ub = inf;
        
    case 13
        % rosenbrock_function
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function 
        func = @rosenbrock_function;

        lb = -inf;
        ub = inf;
        
    case 14
        % powell_function
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @powell_function;

        lb = -inf;
        ub = inf;
        
    case 15
         % quartic_with_noise_function 
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @quartic_with_noise_function;

        lb = -inf;
        ub = inf;
        
    case 16
        % schubert_function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @schubert_function;

        lb = -10;
        ub = 10;
        
    case 17
        % raydan 1 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @raydan_1_function;

        lb = -inf;
        ub = inf;
        
    case 18
        % raydan 2 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @raydan_2_function;

        lb = -inf;
        ub = inf;
        
    case 19
        % extended tridiagonal-1 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_tridiagonal_1_function;

        lb = -inf;
        ub = inf;
        
     case 20
        % extended quadratic penalty qp 1 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_quadratic_penalty_qp1_function;

        lb = -inf;
        ub = inf;
        
    case 21
        % extended quadratic penalty qp 2 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_quadratic_penalty_qp2_function;

        lb = -inf;
        ub = inf;
        
    case 22
        % quadratic qf2 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @quadratic_qf2_function;

        lb = -inf;
        ub = inf;
        
    case 23
        % extended psc1 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_psc1_function;

        lb = -inf;
        ub = inf;
        
    case 24
        % extended bd1 function
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @extended_bd1_function;

        lb = -inf;
        ub = inf;
        
    case 25
        % extended cliff function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_cliff_function;

        lb = -inf;
        ub = inf;
        
    case 26
        % perturbed quadratic diagonal function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @perturbed_quadratic_diagonal_function;

        lb = -inf;
        ub = inf;
%         
    case 27
        % extended hiebert function
        n = 1000;
        x0 = 2000*ones(n,1);

        % Define the objective function
        func = @extended_hiebert_function;

        lb = -inf;
        ub = inf;
        
    case 28
        % extended tet function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_tet_function;

        lb = -inf;
        ub = inf;
        
    case 29
        % diagonal 1 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @diagonal_1_function;

        lb = -inf;
        ub = inf;
%         
    case 30
        % diagonal 3 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @diagonal_3_function;

        lb = -inf;
        ub = inf;
        
    case 31
        % diagonal 5 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @diagonal_5_function;

        lb = -inf;
        ub = inf;
        
    case 32
        % extended maratos function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @extended_maratos_function;

        lb = -inf;
        ub = inf;
        
    case 33
        % eg2 function
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        func = @eg2_function;

        lb = -inf;
        ub = inf;
        
    case 34
        % sinquad function
%         n = 1000;
%         x0 = ones(n,1);
% 
%         % Define the objective function
%         func = @sinquad_function;  
% 
%         lb = -inf;
%         ub = inf;
        
    case 35  
         % griewank_function
        n = 10;
        x0 = 1*ones(n,1);

        % Define the objective function
        func = @griewank_function;

        lb = -inf;
        ub = inf;
        
    case 36  
         % levy13_function
        n = 2;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @levy13_function;

        lb = -inf;
        ub = inf;

    case 37
        % hosaki_function
        n = 2;
        x0 = -2*ones(n,1);

        % Define the objective function
        func = @hosaki_function;

        lb = -inf;
        ub = inf;

    case 38
        % beale_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @beale_function;

        lb = -inf;
        ub = inf;
        
    case 39
        % easom_function
        n = 2;
        x0 = -1*ones(n,1);

        % Define the objective function
        func = @easom_function;

        lb = -inf;
        ub = inf;
        
    case 40
        % price function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @price_function;

        lb = -inf;
        ub = inf;
        
    case 41
        % branin_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @branin_function;

        lb = -inf;
        ub = inf;
        
    case 42
        % trecanni_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @trecanni_function;

        lb = -inf;
        ub = inf;

    case 43
        % booth_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @booth_function;

        lb = -inf;
        ub = inf;
        
    case 44
        % matyas_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @matyas_function;

        lb = -inf;
        ub = inf;
        
    case 45
        % mccormick_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @mccormick_function;

        lb = -3;
        ub = 4;
        
    case 46
        % power_sum_function 
        n = 4;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @power_sum_function;

        lb = -inf;
        ub = inf;
        
    case 47
        % colville_function
        n = 4;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @colville_function;

        lb = -inf;
        ub = inf;
        
    case 48
        % schaffer_function_n_2
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @schaffer_function_n_2;

        lb = -inf;
        ub = inf;

    case 49
        % bohachevsky_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @bohachevsky_function;

        lb = -inf;
        ub = inf;
        
    case 50
        % three_hump_camel_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @three_hump_camel_function;

        lb = -inf;
        ub = inf;
        
    case 51
        % six_hump_camel_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @six_hump_camel_function;

        lb = -inf;
        ub = inf;
        
    case 52
        % drop_wave_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @drop_wave_function;

        lb = -inf;
        ub = inf;
        
    case 53
        % perm_function
        n = 4;
        x0 = ones(n,1);

        % Define the objective function
        func = @perm_function;

        lb = -inf;
        ub = inf;
        
     case 54
        % hartmann_3_d_dimensional_function
        n = 3;
        x0 = ones(n,1);

        % Define the objective function
        func = @hartmann_3_d_dimensional_function;

        lb = -inf;
        ub = inf;
        
     case 55
        % trefethen_4_function
        n = 2;
        x0 = 0.5*ones(n,1);

        % Define the objective function
        func = @trefethen_4_function;

        lb = -inf;
        ub = inf;
        
     case 56
        % zettl_function
        n = 2;
        x0 = 0.1*ones(n,1);

        % Define the objective function
        func = @zettl_function;

        lb = -inf;
        ub = inf;
        
     case 57
        % exp2_function
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @exp2_function;

        lb = -inf;
        ub = inf;
        
     case 58
        % hansen_function
        n = 2;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @hansen_function;

        lb = -inf;
        ub = inf;
        
     case 59
        % schaffer_function_n_4
        n = 2;
        x0 = ones(n,1);

        % Define the objective function
        func = @schaffer_function_n_4;

        lb = -inf;
        ub = inf;
        
     case 60
        % holder_table_function
        n = 2;
        x0 = 10*ones(n,1);

        % Define the objective function
        func = @holder_table_function;

        lb = -inf;
        ub = inf;
        
     case 61
        % gramacy_lee_function
        n = 1;
        x0 = ones(n,1);

        % Define the objective function
        func = @gramacy_lee_function;

        lb = 0.5;
        ub = 3;
        
     case 62
        % eggholder_function
        n = 2;
        x0 = 500*ones(n,1);

        % Define the objective function
        func = @eggholder_function;

        lb = -inf;
        ub = inf;
        
     case 63
        % michalewicz_function 
        n = 2;
        x0 = 2*ones(n,1);

        % Define the objective function
        func = @michalewicz_function;

        lb = -inf;
        ub = inf;
        
    case 64
        % box_betts_exponential_quadratic_sum_function
        n = 3;
        x0 = ones(n,1);

        % Define the objective function
        func = @box_betts_exponential_quadratic_sum_function;

        lb = -inf;
        ub = inf;
        
    case 65
        % cross_in_tray_function
        n=2;
        x0 = ones(n,1);

        % Define the objective function
        func = @cross_in_tray_function;

        lb = -10;
        ub = 10;
        
    case 66
        % himmelblau_function
        n=2;
        x0 = ones(n,1);

        % Define the objective function
        func = @himmelblau_function;

        lb = -inf;
        ub = inf;
        
    case 67
        % forrester_function
        n=1;
        x0 = ones(n,1);

        % Define the objective function
        func = @forrester_function;

        lb = -inf;
        ub = inf;
        
    case 68
        % goldstein-price function
        n = 2;
        x0 = [0;1];

        % Define the objective function
        func = @goldstein_price_function;

        lb = -2;
        ub = 2;
        
end

% test_CNMGE_AG
[f_opt, x_opt,CPU_time, Num_stp] = CNMGE(func, x0, ub, lb);

% Output the global minumum value computed by CNMGE_AG.
fprintf('The global minimum value computed by CNMGE_AG is %12.8f\n',f_opt);

%Output the computational time of CNMGE_AG.
fprintf('The computational time of CNMGE_AG is %8.4f\n',CPU_time);


