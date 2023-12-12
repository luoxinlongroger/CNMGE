clear all
% global test;

% num = [93 107 144 150];

% solveNCPOther = zeros(85,2);
% for test = 145:150
test = input('\n Please input the number of test problem:');

% for i=3:4
%     test = num(i);

switch test
    case 1
        % molecular energy problem
        n = 1000;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 10*ones(n,1);
    case 2
        % ackley_function
        n = 1000;
        xstart = ones(n,1);
        lb = -32.768*ones(n,1);
        ub = 32.768*ones(n,1);

    case 3
        % levy_function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 4
        % schwefel_function
        n = 1000;
        xstart = ones(n,1);
        lb = -500*ones(n,1);
        ub = 500*ones(n,1);

    case 5
        % rastrigin_function
        n = 1000;
        xstart = ones(n,1);
        lb = -5.12*ones(n,1);
        ub = 5.12*ones(n,1);

    case 6
        % styblinski_tang_function
        n = 1000;
        xstart = ones(n,1);
        lb = -5*ones(n,1);
        ub = 5*ones(n,1);

    case 7
        % trid function
        n = 1000;
        xstart = ones(n,1);
        lb = -n^2*ones(n,1);
        ub = n^2*ones(n,1);

    case 8
        % sum_squares_function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 9
        % sphere_function
        n = 1000;
        xstart = ones(n,1);
        lb = -5.12*ones(n,1);
        ub = 5.12*ones(n,1);

    case 10
        % rotated_hyper_ellipsoid_function
        n = 1000;
        xstart = ones(n,1);
        lb = -65.536*ones(n,1);
        ub = 65.536*ones(n,1);

    case 11
        % zakharov_function
        n = 1000;
        xstart = ones(n,1);
        lb = -5*ones(n,1);
        ub = 10*ones(n,1);

    case 12
        % dixon_price_function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 13
        % rosenbrock_function
        n = 1000;
        xstart = 2*ones(n,1);
        lb = -5*ones(n,1);
        ub = 10*ones(n,1);

    case 14
        % powell_function
        n = 1000;
        xstart = ones(n,1);
        lb = -4*ones(n,1);
        ub = 5*ones(n,1);

    case 15
        % quartic_with_noise_function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 16
        % schubert_function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 17
        % raydan 1 function
        n = 1000;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 100*ones(n,1);

    case 18
        % raydan 2 function
        n = 1000;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 100*ones(n,1);

    case 19
        % extended tridiagonal-1 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 20
        % extended quadratic penalty qp 1 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 21
        % extended quadratic penalty qp 2 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 22
        % quadratic qf2 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 23
        % extended psc1 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 24
        % extended bd1 function
        n = 1000;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 25
        % extended cliff function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 26
        % perturbed quadratic diagonal function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 27
        % extended hiebert function -5 5
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 28
        % extended tet function
        n = 1000;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 29
        % diagonal 1 function
        n = 1000;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 30
        % diagonal 3 function
        n = 1000;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 31
        % diagonal 5 function
        n = 1000;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 32
        % extended maratos function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 33
        % eg2 function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 34
        % sinquad function
        n = 1000;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 35
        % griewank_function
        n = 10;
        xstart = ones(n,1);
        lb = -600*ones(n,1);
        ub = 600*ones(n,1);

    case 36
        % levy13_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 37
        % hosaki_function
        n = 2;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 5*ones(n,1);

    case 38
        % beale_function
        n = 2;
        xstart = ones(n,1);
        lb = -4.5*ones(n,1);
        ub = 4.5*ones(n,1);

    case 39
        % easom_function
        n = 2;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 40
        % price function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 41
        % branin_function
        n = 2;
        xstart = ones(n,1);
        lb = -5*ones(n,1);
        ub = 15*ones(n,1);

    case 42
        % trecanni_function
        n = 2;
        xstart = ones(n,1);
        lb = -5*ones(n,1);
        ub = 5*ones(n,1);

    case 43
        % booth_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 44
        % matyas_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 45
        % mccormick_function
        n = 2;
        xstart = ones(n,1);
        lb = -3*ones(n,1);
        ub = 4*ones(n,1);

    case 46
        % power_sum_function
        n = 4;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = n*ones(n,1);

    case 47
        % colville_function
        n = 4;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 48
        % schaffer_function n.2
        n = 2;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 49
        % bohachevsky_function
        n = 2;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 50
        % three_hump_camel_function
        n = 2;
        xstart = ones(n,1);
        lb = -5*ones(n,1);
        ub = 5*ones(n,1);

    case 51
        % six_hump_camel_function
        n = 2;
        xstart = ones(n,1);
        lb = -3*ones(n,1);
        ub = 3*ones(n,1);

    case 52
        % drop_wave_function
        n = 2;
        xstart = ones(n,1);
        lb = -5.12*ones(n,1);
        ub = 5.12*ones(n,1);

    case 53
        % perm_function
        n = 4;
        xstart = ones(n,1);
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 54
        % hartmann_3_d_dimensional_function
        n = 3;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 1*ones(n,1);

    case 55
        % trefethen_4_function
        n = 2;
        xstart = ones(n,1);
        lb = -6.5*ones(n,1);
        ub = 6.5*ones(n,1);

    case 56
        % zettl_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 57
        % exp2_function
        n = 2;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 20*ones(n,1);

    case 58
        % hansen_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 59
        % schaffer_function_n_4
        n = 2;
        xstart = ones(n,1);
        lb = -100*ones(n,1);
        ub = 100*ones(n,1);

    case 60
        % holder_table_function
        n = 2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 61
        % gramacy_lee_function
        n = 1;
        xstart = ones(n,1);
        lb = 0.5*ones(n,1);
        ub = 3*ones(n,1);

    case 62
        % eggholder_function
        n = 2;
        xstart = ones(n,1);
        lb = -512*ones(n,1);
        ub = 512*ones(n,1);

    case 63
        % michalewicz_function
        n = 2;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = pi*ones(n,1);

    case 64
        % box_betts_exponential_quadratic_sum_function
        n = 3;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 65
        % cross_in_tray_function
        n=2;
        xstart = ones(n,1);
        lb = -10*ones(n,1);
        ub = 10*ones(n,1);

    case 66
        % himmelblau_function
        n=2;
        xstart = ones(n,1);
        lb = -6*ones(n,1);
        ub = 6*ones(n,1);

    case 67
        % forrester_function
        n=1;
        xstart = ones(n,1);
        lb = 0*ones(n,1);
        ub = 1*ones(n,1);

    case 68
        % goldstein-price function
        n = 2;
        xstart = ones(n,1);
        lb = -2*ones(n,1);
        ub = 2*ones(n,1);

        %%%%%%%%%%%%%%%%%%%%%%%
    case 69
        n=4;
        xstart = ones(n,1);
        func = 'allinitu';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 71
        n=2;
        xstart = [0;-1];
        func = 'cliff';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 72
        n=10;
        xstart = -ones(n,1);
        func = 'dixon3dq';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 73
        n=2000;
        xstart = zeros(n,1);
        func = 'edensch';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 74
        n=100;
        xstart = zeros(n,1);
        func = 'fletchcr';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 75
        n=500;
        xstart = ones(n,1)/501;
        func = 'genrose';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 76
        n=2;
        xstart = [-5;-7];
        func = 'hairy';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 77
        n=2;
        xstart = [-1.2;1];
        func = 'himmelbb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 78
        n=2;
        xstart = [0.5;0.5];
        func = 'himmelbg';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 79
        n=1000;
        xstart = 1/1001*ones(n,1);
        func = 'indef';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 80
        n=2;
        xstart = [0.3;0.4];
        func = 'jensmp';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 81
        n=1000;
        xstart = 4*ones(n,1);
        func = 'liarwhd';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 82
        n=2;
        xstart = [-500;-700];
        func = 'loghairy';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 83
        n=2;
        xstart = [0;0];
        func = 'maratosb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 84
        n=2;
        xstart = [0.86;0.72];
        func = 'mexhat';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 85
        n=1000;
        xstart = -1*ones(n,1);
        func = 'nondia';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 86
        n=1000;
        xstart = ones(n,1);

        for i=1:n
            if(mod(i,2) == 0)
                xstart = -xstart;
            end
        end
        func = 'nondquar';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 87
        n=1000;
        xstart = 1:n;
        xstart = xstart';

        func = 'penalty1';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);


    case 88
        n=1000;
        xstart = ones(n,1);

        func = 'power1';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);


    case 89
        n=10;
        xstart = ones(n,1);
        func = 'arglinb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 90
        n=10;
        xstart = ones(n,1);
        func = 'arglinc';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 91
        n=1000;
        xstart = ones(n,1);
        func = 'arwhead';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 92
        n=3;
        xstart = ones(n,1);
        func = 'bard';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 93
        n=1000;
        xstart = ones(n,1);
        func = 'bdexp';
        lb = -10*ones(n,1);
        ub = (10+1)*ones(n,1);

    case 94
        n=1000;
        xstart = ones(n,1);
        func = 'bdqrtic';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 95
        n=6;
        xstart = [1 2 1 1 4 3]';
        func = 'biggs6';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 96
        n=3;
        xstart = [0 10 1]';
        func = 'box3';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 97
        n=2;
        xstart = [2 2]';
        func = 'brkmcc';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 98
        n=10;
        xstart = 1/2*ones(n,1);
        func = 'brownal';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 99
        n = 2;
        xstart = ones(n,1);
        func = 'brownbs';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 100
        n = 4;
        xstart = [25 5 -5 -1]';
        func = 'brownden';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 101
        n = 1000;
        xstart = ones(n,1);
        func = 'broydn7d';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 102
        ns = 499;
        n = 2*ns + 2;
        xstart = -2*ones(n,1);
        xstart(1) = -3;
        xstart(2) = -1;
        xstart(3) = -3;
        xstart(4) = -1;
        func = 'chainwoo';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 103
        n = 50;
        xstart = -ones(n,1);
        func = 'chnrosnb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 104
        n = 1000;
        xstart = ones(n,1);
        func = 'cosine';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 105
        n = 1000;
        xstart = 2*ones(n,1);
        xstart(1) = 1;
        func = 'cragglvy';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 106
        n = 2;
        xstart = [-1.2 1]';
        func = 'cube';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 107
        n = 1000;
        xi = 1:n;
        xstart = (0.0001*xi/(n+1))';
        func = 'curly10';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 108
        n = 1000;
        xstart = 3*ones(n,1);
        func = 'dqdrtic';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 109
        n = 1000;
        xstart = 2*ones(n,1);
        func = 'dqrtic';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 110
        n = 1000;
        xstart = zeros(n,1);
        func = 'eg2';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 111
        n = 1000;
        xstart = 2*ones(n,1);
        func = 'engval1';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 112
        n = 3;
        xstart = [1 2 0]';
        func = 'engval2';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 113
        n = 50;
        xstart = -ones(n,1);
        func = 'errinros';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 114
        n = 10;
        xstart = ones(n,1);
        func = 'extrosnb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 115
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        xstart = var/(h);
        func = 'fletcbv3';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 116
        n = 5;
        xstart = 506.2*ones(n,1);
        xstart(1) = -506;
        func = 'genhumps';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 117
        n = 3;
        xstart = [5 2.5 0.15]';
        func = 'gulf';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 118
        n = 50;
        xstart = -3*ones(n,1);
        func = 'hilbertb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 119
        n = 2;
        xstart = [0 2]';
        func = 'himmelbh';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 120
        n = 4;
        xstart = [0.25 0.39 0.415 0.39]';
        func = 'kowosb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 121
        n = 3;
        xstart = [0.02 4000 250]';
        func = 'meyer3';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 122
        n = 1000;
        var = (1:n)';
        xstart = var;
        func = 'noncvxu2';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 123
        n = 5;
        xstart = [0.5 1.5 -1 0.01 0.02]';
        func = 'osbornea';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 124
        n = 100;
        xstart = 1/2*ones(n,1);
        func = 'penalty2';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 125
        n = 1000;
        xstart = 2*ones(n,1);
        func = 'quartc';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 126
        n = 2;
        xstart = [-1.2 1.0]';
        func = 'rosenbr';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 127
        n = 1000;
        var = (1:n)';
        scal = 12.0;
        scale = exp((var-1)*scal/(n-1));
        xstart = 1.0./scale;
        func = 'scosine';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 128
        n = 100;
        sc = 12;
        var1 = (0:n-1)';
        var2 = (1:n)';
        scale = exp(var1*sc/(n-1));
        xstart = 0.0001*(var2.*scale)/(n+1);
        func = 'scurly10';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 129
        n = 2;
        xstart = [4.712389 -1.0]';
        func = 'sineval';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 130
        n = 1000;
        xstart = 0.1*ones(n,1);
        func = 'sinquad';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 131
        n = 2;
        xstart = [1.0 0.1]';
        func = 'sisser';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 132
        n = 1000;
        xstart = ones(n,1);
        xstart(1:2:n) = -1.2*xstart(1:2:n);
        func = 'srosenbr';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 133
        n = 1000;
        xstart = ones(n,1);
        func = 'tridia';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 134
        n = 100;
        var = (1:n)';
        xstart = 1-(var/n);
        func = 'vardim';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 135
        n = 31;
        xstart = zeros(31,1);
        func = 'watson';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 136
        n = 1000;
        xstart = -ones(n,1);
        func = 'woods';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 137
        n = 2;
        xstart = [3 8]';
        func = 'zangwil2';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 138
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        xstart = var*h;
        func = 'fletchbv';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 139
        n = 6;
        xstart = [0 0 1 1 1 1]';
        func = 'heart6ls';
         lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 140
        n = 8;
        xstart = ones(n,1);
        func = 'heart8ls';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 141
        n = 10;
        xstart = ones(n,1);
        xstart(1) = -4;
        xstart(2) = -2;
        func = 'hilberta';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 142
        n = 4;
        xstart = [2.7 90 1500 10]';
        func = 'himmelbf';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 143
        n = 2;
        xstart = [-506.0 -506.2]';
        func = 'humps';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 144
        n = 31;
        xstart = [107.47 0.09203 0.908 102.4 0.1819 ...
            0.8181 97.44 0.284 0.716 96.3 ...
            0.3051 0.6949 93.99 0.3566 0.6434 ...
            89.72 0.468 0.532 83.71 0.6579 ...
            0.3421 78.31 0.8763 0.1237 886.37 ...
            910.01 922.52 926.46 935.56 952.83 ...
            975.73]';
        func = 'methanb8';
         lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 145
        n = 2;
        xstart = [1.0e-30 1.0]';
        func = 'nasty';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 146
        n = 11;
        xstart = [1.3 0.65 0.65 0.7 0.6 ...
            3 5 7 2 4.5 5.5]';
        func = 'osborneb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 147
        n = 3;
        xstart = [0.6 -0.6 20]';
        func = 'yfitu';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 148
        n = 2;
        xstart = [1 1]';
        func = 'denschna';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 149
        n = 2;
        xstart = [1 1]';
        func = 'denschnb';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);

    case 150
        n = 3000;
        xstart = 2*ones(n,1);
        func = 'dixmaank';
        lb = -n*ones(n,1);
        ub = n*ones(n,1);



end
% func = 'test_function';
insigma = [];
inopts.lb = lb;
inopts.ub = ub;
if(n>500)
    inopts.MaxFunEvals = n^2;
end


tic;
[xmin, fmin, counteval, stopflag, out,  bestever ] = cmaes(func,xstart,insigma,inopts);% ,varargin
CMAES_computation_time = toc;

bestever.f

% solveNCPOther(test-68,2) = CMAES_computation_time;
% solveNCPOther(test-68,1) = bestever.f;
% 
% end

% CMAES_result = 'The global minimum value computed by CMA-ES is %12.8f and the computational time of CMA-ES is %6.4f s\n';
% fprintf(CMAES_result,bestever.f,CMAES_computation_time)