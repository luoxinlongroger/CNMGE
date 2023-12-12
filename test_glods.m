clear all;

% solveNCPOther = zeros(85,2);
% for test = 103:133
test = input('\n Please input the number of test problem:');
switch test
    case 1
        % molecular energy problem
        n = 1000;
        x0 = ones(n,1);

        func = @molecular_energy_problem;
        lb = -1e+32;
        ub = 1e+32;

    case 2
        % ackley_function
        n = 1000;
        x0 = ones(n,1);

        func = @ackley_function;
        lb = -1e+32;
        ub = 1e+32;

    case 3
        % levy_function
        n = 1000;
        x0 = 2*ones(n,1);

        func = @levy_function;
        lb = -1e+32;
        ub = 1e+32;

    case 4
        % schwefel_function
        n = 1000;
        x0 = 200*ones(n,1);

        func = @schwefel_function;
        lb = -500;
        ub = 500;

    case 5
        % rastrigin_function
        n = 1000;
        x0 = ones(n,1);

        func = @rastrigin_function;
        lb = -1e+32;
        ub = 1e+32;

    case 6
        % styblinski_tang_function
        n = 1000;
        x0 = ones(n,1);

        func = @styblinski_tang_function;
        lb = -1e+32;
        ub = 1e+32;

    case 7
        % trid function
        n = 1000;
        x0 = ones(n,1);

        func = @trid_function;
        lb = -1e+32;
        ub = 1e+32;

    case 8
        % sum_squares_function
        n = 1000;
        x0 = ones(n,1);

        func = @sum_squares_function;
        lb = -1e+32;
        ub = 1e+32;

    case 9
        % sphere_function
        n = 1000;
        x0 = ones(n,1);

        func = @sphere_function;
        lb = -1e+32;
        ub = 1e+32;

    case 10
        % rotated_hyper_ellipsoid_function
        n = 1000;
        x0 = ones(n,1);

        func = @rotated_hyper_ellipsoid_function;
        lb = -1e+32;
        ub = 1e+32;

    case 11
        % zakharov_function
        n = 1000;
        x0 = (5e-5)*ones(n,1);

        func = @zakharov_function;
        lb = -1e+32;
        ub = 1e+32;

    case 12
        % dixon_price_function
        n = 1000;
        x0 = ones(n,1);

        func = @dixon_price_function;
        lb = -1e+32;
        ub = 1e+32;

    case 13
        % rosenbrock_function
        n = 1000;
        x0 = 2*ones(n,1);

        func = @rosenbrock_function;
        lb = -1e+32;
        ub = 1e+32;

    case 14
        % powell_function
        n = 1000;
        x0 = 2*ones(n,1);

        func = @powell_function;
        lb = -1e+32;
        ub = 1e+32;

    case 15
        % quartic_with_noise_function
        n = 1000;
        x0 = ones(n,1);

        func = @quartic_with_noise_function;
        lb = -1e+32;
        ub = 1e+32;

    case 16
        % schubert_function
        n = 1000;
        x0 = ones(n,1);

        func = @schubert_function;

        lb = -10;
        ub = 10;

    case 17
        % raydan 1 function
        n = 1000;
        x0 = ones(n,1);

        func = @raydan_1_function;
        lb = -1e+32;
        ub = 1e+32;

    case 18
        % raydan 2 function
        n = 1000;
        x0 = ones(n,1);

        func = @raydan_2_function;
        lb = -1e+32;
        ub = 1e+32;

    case 19
        % extended tridiagonal-1 function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_tridiagonal_1_function;
        lb = -1e+32;
        ub = 1e+32;

    case 20
        % extended quadratic penalty qp 1 function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_quadratic_penalty_qp1_function;
        lb = -1e+32;
        ub = 1e+32;

    case 21
        % extended quadratic penalty qp 2 function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_quadratic_penalty_qp2_function;
        lb = -1e+32;
        ub = 1e+32;

    case 22
        % quadratic qf2 function
        n = 1000;
        x0 = ones(n,1);

        func = @quadratic_qf2_function;
        lb = -1e+32;
        ub = 1e+32;

    case 23
        % extended psc1 function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_psc1_function;
        lb = -1e+32;
        ub = 1e+32;

    case 24
        % extended bd1 function
        n = 1000;
        x0 = 2*ones(n,1);

        func = @extended_bd1_function;
        lb = -1e+32;
        ub = 1e+32;

    case 25
        % extended cliff function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_cliff_function;
        lb = -1e+32;
        ub = 1e+32;

    case 26
        % perturbed quadratic diagonal function
        n = 1000;
        x0 = ones(n,1);

        func = @perturbed_quadratic_diagonal_function;
        lb = -1e+32;
        ub = 1e+32;

    case 27
        % extended hiebert function
        n = 1000;
        x0 = 2000*ones(n,1);

        func = @extended_hiebert_function;
        lb = -1e+32;
        ub = 1e+32;

    case 28
        % extended tet function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_tet_function;
        lb = -1e+32;
        ub = 1e+32;

    case 29
        % diagonal 1 function
        n = 1000;
        x0 = ones(n,1);

        func = @diagonal_1_function;

        lb = -1e+32;
        ub = 1e+32;

    case 30
        % diagonal 3 function
        n = 1000;
        x0 = ones(n,1);

        func = @diagonal_3_function;
        lb = -1e+32;
        ub = 1e+32;

    case 31
        % diagonal 5 function
        n = 1000;
        x0 = ones(n,1);

        func = @diagonal_5_function;
        lb = -1e+32;
        ub = 1e+32;

    case 32
        % extended maratos function
        n = 1000;
        x0 = ones(n,1);

        func = @extended_maratos_function;
        lb = -1e+32;
        ub = 1e+32;

    case 33
        % eg2 function
        n = 1000;
        x0 = ones(n,1);

        func = @eg2_function;
        lb = -1e+32;
        ub = 1e+32;

    case 34
        % sinquad function
        n = 1000;
        x0 = ones(n,1);

        func = @sinquad_function;
        lb = -1e+32;
        ub = 1e+32;

    case 35
        % griewank_function
        n = 10;
        x0 = dlarray(1*ones(n,1));

        func = @griewank_function;
        lb = -1e+32;
        ub = 1e+32;

    case 36
        % levy13_function
        n = 2;
        x0 = dlarray(2*ones(n,1));

        func = @levy13_function;
        lb = -1e+32;
        ub = 1e+32;

    case 37
        % hosaki_function
        n = 2;
        x0 = dlarray(-2*ones(n,1));

        func = @hosaki_function;
        lb = 0;
        ub = 5;

    case 38
        % beale_function
        n = 2;
        x0 = ones(n,1);

        func = @beale_function;
        lb = -1e+32;
        ub = 1e+32;

    case 39
        % easom_function
        n = 2;
        x0 = -1*ones(n,1);

        func = @easom_function;
        lb = -1e+32;
        ub = 1e+32;

    case 40
        % price function
        n = 2;
        x0 = ones(n,1);

        func = @price_function;
        lb = -1e+32;
        ub = 1e+32;

    case 41
        % branin_function
        n = 2;
        x0 = ones(n,1);

        func = @branin_function;
        lb = -1e+32;
        ub = 1e+32;

    case 42
        % trecanni_function
        n = 2;
        x0 = ones(n,1);

        func = @trecanni_function;
        lb = -1e+32;
        ub = 1e+32;

    case 43
        % booth_function
        n = 2;
        x0 = ones(n,1);

        func = @booth_function;
        lb = -1e+32;
        ub = 1e+32;

    case 44
        % matyas_function
        n = 2;
        x0 = ones(n,1);

        func = @matyas_function;
        lb = -1e+32;
        ub = 1e+32;

    case 45
        % mccormick_function
        n = 2;
        x0 = ones(n,1);

        func = @mccormick_function;
        lb = -3;
        ub = 4;

    case 46
        % power_sum_function
        n = 4;
        x0 = 2*ones(n,1);

        func = @power_sum_function;
        lb = -1e+32;
        ub = 1e+32;

    case 47
        % colville_function
        n = 4;
        x0 = 2*ones(n,1);

        func = @colville_function;
        lb = -1e+32;
        ub = 1e+32;

    case 48
        % schaffer_function_n_2
        n = 2;
        x0 = ones(n,1);

        func = @schaffer_function_n_2;
        lb = -1e+32;
        ub = 1e+32;

    case 49
        % bohachevsky_function
        n = 2;
        x0 = ones(n,1);

        func = @bohachevsky_function;
        lb = -1e+32;
        ub = 1e+32;

    case 50
        % three_hump_camel_function
        n = 2;
        x0 = ones(n,1);

        func = @three_hump_camel_function;
        lb = -1e+32;
        ub = 1e+32;

    case 51
        % six_hump_camel_function
        n = 2;
        x0 = ones(n,1);

        func = @six_hump_camel_function;
        lb = -1e+32;
        ub = 1e+32;

    case 52
        % drop_wave_function
        n = 2;
        x0 = ones(n,1);

        func = @drop_wave_function;
        lb = -1e+32;
        ub = 1e+32;

    case 53
        % perm_function
        n = 4;
        x0 = ones(n,1);

        func = @perm_function;
        lb = -1e+32;
        ub = 1e+32;

    case 54
        % hartmann_3_d_dimensional_function
        n = 3;
        x0 = ones(n,1);

        func = @hartmann_3_d_dimensional_function;
        lb = -1e+32;
        ub = 1e+32;

    case 55
        % trefethen_4_function
        n = 2;
        x0 = 0.5*ones(n,1);

        func = @trefethen_4_function;
        lb = -1e+32;
        ub = 1e+32;

    case 56
        % zettl_function
        n = 2;
        x0 = 0.1*ones(n,1);

        func = @zettl_function;
        lb = -1e+32;
        ub = 1e+32;

    case 57
        % exp2_function
        n = 2;
        x0 = ones(n,1);

        func = @exp2_function;
        lb = -1e+32;
        ub = 1e+32;

    case 58
        % hansen_function
        n = 2;
        x0 = 2*ones(n,1);

        func = @hansen_function;
        lb = -1e+32;
        ub = 1e+32;

    case 59
        % schaffer_function_n_4
        n = 2;
        x0 = ones(n,1);
        func = @schaffer_function_n_4;

        lb = -1e+32;
        ub = 1e+32;

    case 60
        % holder_table_function
        n = 2;
        x0 = 10*ones(n,1);
        func = @holder_table_function;

        lb = -10;
        ub = 10;

    case 61
        % gramacy_lee_function
        n = 1;
        x0 = ones(n,1);

        func = @gramacy_lee_function;

        lb = 0.5;
        ub = 3;

    case 62
        % eggholder_function
        n = 2;
        x0 = 500*ones(n,1);

        func = @eggholder_function;

        lb = -512;
        ub = 512;

    case 63
        % michalewicz_function
        n = 2;
        x0 = 2*ones(n,1);

        func = @michalewicz_function;
        lb = -1e+32;
        ub = 1e+32;

    case 64
        % box_betts_exponential_quadratic_sum_function
        n = 3;
        x0 = ones(n,1);
        func = @box_betts_exponential_quadratic_sum_function;

        lb = -1e+32;
        ub = 1e+32;

    case 65
        % cross_in_tray_function
        n=2;
        x0 = ones(n,1);

        func = @cross_in_tray_function;

        lb = -10;
        ub = 10;


    case 66
        % himmelblau_function
        n=2;
        x0 = ones(n,1);

        func = @himmelblau_function;

        lb = -1e+32;
        ub = 1e+32;

    case 67
        % forrester_function
        n=1;
        x0 = ones(n,1);

        func = @forrester_function;
        lb = 0;
        ub = 1;

    case 68
        % goldstein-price function
        n = 2;
        x0 = [0;1];
        func = @goldstein_price_function;

        lb = -2;
        ub = 2;

    case 69
        n=4;
        x0 = ones(n,1);
        func = @allinitu;
        lb = -1e+32;
        ub = 1e+32;

    case 71
        n=2;
        x0 = [0;-1];
        func = @cliff;
        lb = -1e+32;
        ub = 1e+32;

    case 72
        n=10;
        x0 = -ones(n,1);
        func = @dixon3dq;
        lb = -1e+32;
        ub = 1e+32;

    case 73
        n=2000;
        x0 = zeros(n,1);
        func = @edensch;
        lb = -1e+32;
        ub = 1e+32;

    case 74
        n=100;
        x0 = zeros(n,1);
        func = @fletchcr;
        lb = -1e+32;
        ub = 1e+32;

    case 75
        n=500;
        x0 = ones(n,1)/501;
        func = @genrose;
        lb = -1e+32;
        ub = 1e+32;

    case 76
        n=2;
        x0 = [-5;-7];
        func = @hairy;
        lb = -1e+32;
        ub = 1e+32;

    case 77         %卡住算不了
        n=2;
        x0 = [-1.2;1];
        func = @himmelbb;
        lb = -1e+32;
        ub = 1e+32;

    case 78
        n=2;
        x0 = [0.5;0.5];
        func = @himmelbg;
        lb = -1e+32;
        ub = 1e+32;

    case 79
        n=1000;
        x0 = 1/1001*ones(n,1);
        func = @indef;
        lb = -1e+32;
        ub = 1e+32;

    case 80
        n=2;
        x0 = [0.3;0.4];
        func = @jensmp;
        lb = -1e+32;
        ub = 1e+32;

    case 81
        n=1000;
        x0 = 4*ones(n,1);
        func = @liarwhd;
        lb = -1e+32;
        ub = 1e+32;

    case 82
        n=2;
        x0 = [-500;-700];
        func = @loghairy;
        lb = -1e+32;
        ub = 1e+32;

    case 83
        n=2;
        x0 = [0;0];
        func = @maratosb;
        lb = -1e+32;
        ub = 1e+32;

    case 84
        n=2;
        x0 = [0.86;0.72];
        func = @mexhat;
        lb = -1e+32;
        ub = 1e+32;

    case 85
        n=1000;
        x0 = -1*ones(n,1);
        func = @nondia;
        lb = -1e+32;
        ub = 1e+32;

    case 86
        n=1000;
        x0 = ones(n,1);

        for i=1:n
            if(mod(i,2) == 0)
                x0 = -x0;
            end
        end
        func = @nondquar;
        lb = -1e+32;
        ub = 1e+32;

    case 87
        n=1000;
        x0 = 1:n;
        x0 = x0';

        func = @penalty1;
        lb = -1e+32;
        ub = 1e+32;


    case 88
        n=1000;
        x0 = ones(n,1);

        func = @power1;
        lb = -1e+32;
        ub = 1e+32;


    case 89
        n=10;
        x0 = ones(n,1);
        func = @arglinb;
        lb = -1e+32;
        ub = 1e+32;

    case 90
        n=10;
        x0 = ones(n,1);
        func = @arglinc;
        lb = -1e+32;
        ub = 1e+32;

    case 91
        n=1000;
        x0 = ones(n,1);
        func = @arwhead;
        lb = -1e+32;
        ub = 1e+32;

    case 92
        n=3;
        x0 = ones(n,1);
        func = @bard;
        lb = -1e+32;
        ub = 1e+32;

    case 93
        n=1000;
        x0 = ones(n,1);
        func = @bdexp;
        lb = -1e+32;
        ub = 1e+32;

    case 94
        n=1000;
        x0 = ones(n,1);
        func = @bdqrtic;
        lb = -1e+32;
        ub = 1e+32;

    case 95
        n=6;
        x0 = [1 2 1 1 4 3]';
        func = @biggs6;
        lb = -1e+32;
        ub = 1e+32;

    case 96
        n=3;
        x0 = [0 10 1]';
        func = @box3;
        lb = -1e+32;
        ub = 1e+32;

    case 97
        n=2;
        x0 = [2 2]';
        func = @brkmcc;
        lb = -1e+32;
        ub = 1e+32;

    case 98
        n=10;
        x0 = 1/2*ones(n,1);
        func = @brownal;
        lb = -1e+32;
        ub = 1e+32;

    case 99
        n = 2;
        x0 = ones(n,1);
        func = @brownbs;
        lb = -1e+32;
        ub = 1e+32;

    case 100
        n = 4;
        x0 = [25 5 -5 -1]';
        func = @brownden;
        lb = -1e+32;
        ub = 1e+32;

    case 101
        n = 1000;
        x0 = ones(n,1);
        func = @broydn7d;
        lb = -1e+32;
        ub = 1e+32;

    case 102
        ns = 499;
        n = 2*ns + 2;
        x0 = -2*ones(n,1);
        x0(1) = -3;
        x0(2) = -1;
        x0(3) = -3;
        x0(4) = -1;
        func = @chainwoo;
        lb = -1e+32;
        ub = 1e+32;

    case 103
        n = 50;
        x0 = -ones(n,1);
        func = @chnrosnb;
        lb = -1e+32;
        ub = 1e+32;

    case 104
        n = 1000;
        x0 = ones(n,1);
        func = @cosine;
        lb = -1e+32;
        ub = 1e+32;

    case 105
        n = 1000;
        x0 = 2*ones(n,1);
        x0(1) = 1;
        func = @cragglvy;
        lb = -1e+32;
        ub = 1e+32;

    case 106
        n = 2;
        x0 = [-1.2 1]';
        func = @cube;
        lb = -1e+32;
        ub = 1e+32;

    case 107
        n = 1000;
        xi = 1:n;
        x0 = (0.0001*xi/(n+1))';
        func = @curly10;
        lb = -1e+32;
        ub = 1e+32;

    case 108
        n = 1000;
        x0 = 3*ones(n,1);
        func = @dqdrtic;
        lb = -1e+32;
        ub = 1e+32;

    case 109
        n = 1000;
        x0 = 2*ones(n,1);
        func = @dqrtic;
        lb = -1e+32;
        ub = 1e+32;

    case 110
        n = 1000;
        x0 = zeros(n,1);
        func = @eg2;
        lb = -1e+32;
        ub = 1e+32;

    case 111
        n = 1000;
        x0 = 2*ones(n,1);
        func = @engval1;
        lb = -1e+32;
        ub = 1e+32;

    case 112
        n = 3;
        x0 = [1 2 0]';
        func = @engval2;
        lb = -1e+32;
        ub = 1e+32;

    case 113
        n = 50;
        x0 = -ones(n,1);
        func = @errinros;
        lb = -1e+32;
        ub = 1e+32;

    case 114
        n = 10;
        x0 = ones(n,1);
        func = @extrosnb;
        lb = -1e+32;
        ub = 1e+32;

    case 115
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var/(h);
        func = @fletcbv3;
        lb = -1e+32;
        ub = 1e+32;

    case 116
        n = 5;
        x0 = 506.2*ones(n,1);
        x0(1) = -506;
        func = @genhumps;
        lb = -1e+32;
        ub = 1e+32;

    case 117
        n = 3;
        x0 = [5 2.5 0.15]';
        func = @gulf;
        lb = -1e+32;
        ub = 1e+32;

    case 118
        n = 50;
        x0 = -3*ones(n,1);
        func = @hilbertb;
        lb = -1e+32;
        ub = 1e+32;

    case 119
        n = 2;
        x0 = [0 2]';
        func = @himmelbh;
        lb = -1e+32;
        ub = 1e+32;

    case 120
        n = 4;
        x0 = [0.25 0.39 0.415 0.39]';
        func = @kowosb;
        lb = -1e+32;
        ub = 1e+32;

    case 121
        n = 3;
        x0 = [0.02 4000 250]';
        func = @meyer3;
        lb = -1e+32;
        ub = 1e+32;

    case 122
        n = 1000;
        var = (1:n)';
        x0 = var;
        func = @noncvxu2;
        lb = -1e+32;
        ub = 1e+32;

    case 123
        n = 5;
        x0 = [0.5 1.5 -1 0.01 0.02]';
        func = @osbornea;
        lb = -1e+32;
        ub = 1e+32;

    case 124
        n = 100;
        x0 = 1/2*ones(n,1);
        func = @penalty2;
        lb = -1e+32;
        ub = 1e+32;

    case 125
        n = 1000;
        x0 = 2*ones(n,1);
        func = @quartc;
        lb = -1e+32;
        ub = 1e+32;

    case 126
        n = 2;
        x0 = [-1.2 1.0]';
        func = @rosenbr;
        lb = -1e+32;
        ub = 1e+32;

    case 127
        n = 1000;
        var = (1:n)';
        scal = 12.0;
        scale = exp((var-1)*scal/(n-1));
        x0 = 1.0./scale;
        func = @scosine;
        lb = -1e+32;
        ub = 1e+32;

    case 128
        n = 100;
        sc = 12;
        var1 = (0:n-1)';
        var2 = (1:n)';
        scale = exp(var1*sc/(n-1));
        x0 = 0.0001*(var2.*scale)/(n+1);
        func = @scurly10;
        lb = -1e+32;
        ub = 1e+32;

    case 129
        n = 2;
        x0 = [4.712389 -1.0]';
        func = @sineval;
        lb = -1e+32;
        ub = 1e+32;

    case 130
        n = 1000;
        x0 = 0.1*ones(n,1);
        func = @sinquad;
        lb = -1e+32;
        ub = 1e+32;

    case 131
        n = 2;
        x0 = [1.0 0.1]';
        func = @sisser;
        lb = -1e+32;
        ub = 1e+32;

    case 132
        n = 1000;
        x0 = ones(n,1);
        x0(1:2:n) = -1.2*x0(1:2:n);
        func = @srosenbr;
        lb = -1e+32;
        ub = 1e+32;

    case 133
        n = 1000;
        x0 = ones(n,1);
        func = @tridia;
        lb = -1e+32;
        ub = 1e+32;

    case 134
        n = 100;
        var = (1:n)';
        x0 = 1-(var/n);
        func = @vardim;
        lb = -1e+32;
        ub = 1e+32;

    case 135
        n = 31;
        x0 = zeros(31,1);
        func = @watson;
        lb = -1e+32;
        ub = 1e+32;

    case 136
        n = 1000;
        x0 = -ones(n,1);
        func = @woods;
        lb = -1e+32;
        ub = 1e+32;

    case 137
        n = 2;
        x0 = [3 8]';
        func = @zangwil2;
        lb = -1e+32;
        ub = 1e+32;

    case 138
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var*h;
        func = @fletchbv;
        lb = -1e+32;
        ub = 1e+32;

    case 139
        n = 6;
        x0 = [0 0 1 1 1 1]';
        func = @heart6ls;
        lb = -1e+32;
        ub = 1e+32;

    case 140
        n = 8;
        x0 = ones(n,1);
        func = @heart8ls;
        lb = -1e+32;
        ub = 1e+32;

    case 141
        n = 10;
        x0 = ones(n,1);
        x0(1) = -4;
        x0(2) = -2;
        func = @hilberta;
        lb = -1e+32;
        ub = 1e+32;

    case 142
        n = 4;
        x0 = [2.7 90 1500 10]';
        func = @himmelbf;
        lb = -1e+32;
        ub = 1e+32;

    case 143
        n = 2;
        x0 = [-506.0 -506.2]';
        func = @humps;
        lb = -1e+32;
        ub = 1e+32;

    case 144
        n = 31;
        x0 = [107.47 0.09203 0.908 102.4 0.1819 ...
            0.8181 97.44 0.284 0.716 96.3 ...
            0.3051 0.6949 93.99 0.3566 0.6434 ...
            89.72 0.468 0.532 83.71 0.6579 ...
            0.3421 78.31 0.8763 0.1237 886.37 ...
            910.01 922.52 926.46 935.56 952.83 ...
            975.73]';
        func = @methanb8;
        lb = -1e+32;
        ub = 1e+32;

    case 145
        n = 2;
        x0 = [1.0e-30 1.0]';
        func = @nasty;
        lb = -1e+32;
        ub = 1e+32;

    case 146
        n = 11;
        x0 = [1.3 0.65 0.65 0.7 0.6 ...
            3 5 7 2 4.5 5.5]';
        func = @osborneb;
        lb = -1e+32;
        ub = 1e+32;

    case 147
        n = 3;
        x0 = [0.6 -0.6 20]';
        func = @yfitu;
        lb = -1e+32;
        ub = 1e+32;

    case 148
        n = 2;
        x0 = [1 1]';
        func = @denschna;
        lb = -1e+32;
        ub = 1e+32;

    case 149
        n = 2;
        x0 = [1 1]';
        func = @denschnb;
        lb = -1e+32;
        ub = 1e+32;
        

    case 150
        n = 3000;
        x0 = 2*ones(n,1);
        func = @dixmaank;
        lb = -1e+32;
        ub = 1e+32;

end



lbound = lb*ones(n,1);
ubound = ub*ones(n,1);

tic
[glods_profile,Plist,flist,alfa,radius,func_eval,outmin] = glods(func,[],x0,lbound,ubound);
time = toc;


% solveNCPOther(test-68,2) = time;
% solveNCPOther(test-68,1) = outmin;
% 
% end