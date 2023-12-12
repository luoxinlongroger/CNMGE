clear; 

test = input('\n Please input the number of test problem:');
switch test
    case 69
        n=4;
        x0 = ones(n,1);
        func = @allinitu;

    case 71
        n=2;
        x0 = [0;-1];
        func = @cliff;

    case 72
        n=10;
        x0 = -ones(n,1);
        func = @dixon3dq;

    case 73
        n=2000;
        x0 = zeros(n,1);
        func = @edensch;

    case 74
        n=100;
        x0 = zeros(n,1);
        func = @fletchcr;

    case 75
        n=500;
        x0 = ones(n,1)/501;
        func = @genrose;

    case 76
        n=2;
        x0 = [-5;-7];
        func = @hairy;

    case 77
        n=2;
        x0 = [-1.2;1];
        func = @himmelbb;

    case 78
        n=2;
        x0 = [0.5;0.5];
        func = @himmelbg;

    case 79
        n=1000;
        x0 = 1/1001*ones(n,1);
        func = @indef;

    case 80
        n=2;
        x0 = [0.3;0.4];
        func = @jensmp;

    case 81
        n=1000;
        x0 = 4*ones(n,1);
        func = @liarwhd;

    case 82
        n=2;
        x0 = [-500;-700];
        func = @loghairy;

    case 83
        n=2;
        x0 = [0;0];
        func = @maratosb;

    case 84
        n=2;
        x0 = [0.86;0.72];
        func = @mexhat;

    case 85
        n=1000;
        x0 = -1*ones(n,1);
        func = @nondia;

    case 86
        n=1000;
        x0 = ones(n,1);

        for i=1:n
            if(mod(i,2) == 0)
                x0 = -x0;
            end
        end
        func = @nondquar;

    case 87
        n=1000;
        x0 = 1:n;
        x0 = x0';

        func = @penalty1;


    case 88
        n=1000;
        x0 = ones(n,1);

        func = @power1;


    case 89
        n=10;
        x0 = ones(n,1);
        func = @arglinb;

    case 90 
        n=10;
        x0 = ones(n,1);
        func = @arglinc;   

    case 91
        n=1000;
        x0 = ones(n,1);
        func = @arwhead;

    case 92
        n=3;
        x0 = ones(n,1);
        func = @bard;

    case 93
        n=1000;
        x0 = ones(n,1);
        func = @bdexp;

    case 94
        n=1000;
        x0 = ones(n,1);
        func = @bdqrtic;

    case 95
        n=6;
        x0 = [1 2 1 1 4 3]';
        func = @biggs6;

    case 96
        n=3;
        x0 = [0 10 1]';
        func = @box3;

    case 97
        n=2;
        x0 = [2 2]';
        func = @brkmcc;

    case 98
        n=10;
        x0 = 1/2*ones(n,1);
        func = @brownal;

    case 99
        n = 2; 
        x0 = ones(n,1);
        func = @brownbs;

    case 100
        n = 4;
        x0 = [25 5 -5 -1]';
        func = @brownden;

    case 101
        n = 1000;
        x0 = ones(n,1);
        func = @broydn7d;

    case 102
        ns = 499;
        n = 2*ns + 2;
        x0 = -2*ones(n,1);
        x0(1) = -3;
        x0(2) = -1;
        x0(3) = -3;
        x0(4) = -1;
        func = @chainwoo;

    case 103
        n = 50;
        x0 = -ones(n,1);
        func = @chnrosnb;

    case 104
        n = 1000;
        x0 = ones(n,1);
        func = @cosine;

    case 105
        n = 1000;
        x0 = 2*ones(n,1);
        x0(1) = 1;
        func = @cragglvy;

    case 106
        n = 2;
        x0 = [-1.2 1]';
        func = @cube;

    case 107
        n = 1000;
        xi = 1:n;
        x0 = (0.0001*xi/(n+1))';
        func = @curly10;

    case 108
        n = 1000;
        x0 = 3*ones(n,1);
        func = @dqdrtic;

    case 109
        n = 1000;
        x0 = 2*ones(n,1);
        func = @dqrtic;

    case 110
        n = 1000;
        x0 = zeros(n,1);
        func = @eg2;

    case 111
        n = 1000;
        x0 = 2*ones(n,1);
        func = @engval1;

    case 112
        n = 3;
        x0 = [1 2 0]';
        func = @engval2;

    case 113
        n = 50;
        x0 = -ones(n,1);
        func = @errinros;

    case 114
        n = 10;
        x0 = ones(n,1);
        func = @extrosnb;

    case 115
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var/(h);
        func = @fletcbv3;

    case 116
        n = 5;
        x0 = 506.2*ones(n,1);
        x0(1) = -506;
        func = @genhumps;

    case 117
        n = 3;
        x0 = [5 2.5 0.15]';
        func = @gulf;

    case 118
        n = 50;
        x0 = -3*ones(n,1);
        func = @hilbertb;

    case 119
        n = 2;
        x0 = [0 2]';
        func = @himmelbh;

    case 120
        n = 4;
        x0 = [0.25 0.39 0.415 0.39]';
        func = @kowosb;

    case 121
        n = 3;
        x0 = [0.02 4000 250]';
        func = @meyer3;

    case 122
        n = 1000;
        var = (1:n)';
        x0 = var;
        func = @noncvxu2;

    case 123
        n = 5;
        x0 = [0.5 1.5 -1 0.01 0.02]';
        func = @osbornea;

    case 124
        n = 100;
        x0 = 1/2*ones(n,1);
        func = @penalty2;

    case 125
        n = 1000;
        x0 = 2*ones(n,1);
        func = @quartc;

    case 126
        n = 2;
        x0 = [-1.2 1.0]';
        func = @rosenbr;

    case 127
        n = 1000;
        var = (1:n)';
        scal = 12.0;
        scale = exp((var-1)*scal/(n-1));
        x0 = 1.0./scale;
        func = @scosine;

    case 128
        n = 100;
        sc = 12;
        var1 = (0:n-1)';
        var2 = (1:n)';
        scale = exp(var1*sc/(n-1));
        x0 = 0.0001*(var2.*scale)/(n+1);
        func = @scurly10;     

    case 129
        n = 2;
        x0 = [4.712389 -1.0]';
        func = @sineval;  

    case 130
        n = 1000;
        x0 = 0.1*ones(n,1);
        func = @sinquad;  

    case 131
        n = 2;
        x0 = [1.0 0.1]';
        func = @sisser;  

    case 132
        n = 1000;
        x0 = ones(n,1);
        x0(1:2:n) = -1.2*x0(1:2:n);
        func = @srosenbr;  

    case 133
        n = 1000;
        x0 = ones(n,1);
        func = @tridia;  

    case 134
        n = 100;
        var = (1:n)';
        x0 = 1-(var/n);
        func = @vardim;  

    case 135
        n = 31;
        x0 = zeros(31,1);
        func = @watson;  

    case 136
        n = 1000;
        x0 = -ones(n,1);
        func = @woods;  

    case 137
        n = 2;
        x0 = [3 8]';
        func = @zangwil2; 

    case 138
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var*h;
        func = @fletchbv;

    case 139
        n = 6;
        x0 = [0 0 1 1 1 1]';
        func = @heart6ls;

    case 140
        n = 8;
        x0 = ones(n,1);
        func = @heart8ls;

    case 141
        n = 10;
        x0 = ones(n,1);
        x0(1) = -4;
        x0(2) = -2;
        func = @hilberta;

    case 142
        n = 4;
        x0 = [2.7 90 1500 10]';
        func = @himmelbf;

    case 143
        n = 2;
        x0 = [-506.0 -506.2]';
        func = @humps;

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

    case 145
        n = 2;
        x0 = [1.0e-30 1.0]';
        func = @nasty;

    case 146
        n = 11;
        x0 = [1.3 0.65 0.65 0.7 0.6 ...
              3 5 7 2 4.5 5.5]';
        func = @osborneb;

    case 147
        n = 3;
        x0 = [0.6 -0.6 20]';
        func = @yfitu;

    case 148
        n = 2;
        x0 = [1 1]';
        func = @denschna;

    case 149
        n = 2;
        x0 = [1 1]';
        func = @denschnb;

    case 150
        n = 3000;
        x0 = 2*ones(n,1);
        func = @dixmaank;
end

tic
[x,f,info] = VRBBO(func,x0,[], []);

toc
f
