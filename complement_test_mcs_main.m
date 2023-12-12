%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% runmcs.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test driver for running MCS 
% with advice on how to choose the tuning parameters
%
% applied to the Jones test functions with default boxes
%
% To solve your own problem, copy this file (to ownmcs.m, say)
% and adapt the first part labelled `problem definition'.
% Then run the driver by typing `ownmcs' at the Matlab prompt.
%
% If you are not satisfied with the results, or if you want to play 
% with the tuning parameters, you also need to adapt the second part
% labelled `change MCS settings'. In particular, for wide bounds,
% you'll probably get better results if you supply your own 
% initialization list.
% 
% On typing `runmcs' at the Matlab prompt,
% the unmodified file produces test results for the six-hump camel
% function; by only changing the data assignment you can also get
% results for the other test functions from Jones et al.
% You may also play with the bounds by modifying the default bounds.
% 

clear; clear mex; clear global; 
format compact;format long

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% problem definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define objective function
%
% Each function value f=f(x) is obtained by MCS using a call 
% f=feval(fcn,data,x)
% where data is an arbitrary data object passed to mcs. If you have 
% several data items, put them into a cell array or structure array
% called data.
% If there are no data, use fcn='feval' and the function name as data.
%
fcn = 'feval';

% num = [75 94 101 107 108 109 110 111 115 127 130 132 133 136 138 150];
% num = [115 127 130 132 133 136 138 150];
num = [93  107];

solveNCPOther = zeros(85,2);
% for test = 92:92
% for i=1:8
%     test = num(i);

test = input('\n Please input the number of test problem:');
switch test
    
    case 69
        n=4;
        x0 = ones(n,1);
        data = 'allinitu';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 71
        n=2;
        x0 = [0;-1];
        data = 'cliff';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 72
        n=10;
        x0 = -ones(n,1);
        data = 'dixon3dq';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 73
        n=2000;
        x0 = zeros(n,1);
        data = 'edensch';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 74
        n=100;
        x0 = zeros(n,1);
        data = 'fletchcr';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 75
        n=500;
        x0 = ones(n,1)/501;
        data = 'genrose';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 76
        n=2;
        x0 = [-5;-7];
        data = 'hairy';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 77         %卡住算不了
        n=2;
        x0 = [-1.2;1];
        data = 'himmelbb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 78
        n=2;
        x0 = [0.5;0.5];
        data = 'himmelbg';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 79
        n=1000;
        x0 = 1/1001*ones(n,1);
        data = 'indef';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 80
        n=2;
        x0 = [0.3;0.4];
        data = 'jensmp';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 81
        n=1000;
        x0 = 4*ones(n,1);
        data = 'liarwhd';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 82
        n=2;
        x0 = [-500;-700];
        data = 'loghairy';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 83
        n=2;
        x0 = [0;0];
        data = 'maratosb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 84
        n=2;
        x0 = [0.86;0.72];
        data = 'mexhat';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 85
        n=1000;
        x0 = -1*ones(n,1);
        data = 'nondia';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 86
        n=1000;
        x0 = ones(n,1);

        for i=1:n
            if(mod(i,2) == 0)
                x0 = -x0;
            end
        end
        data = 'nondquar';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 87
        n=1000;
        x0 = 1:n;
        x0 = x0';

        data = 'penalty1';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);


    case 88
        n=1000;
        x0 = ones(n,1);

        data = 'power1';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);


    case 89
        n=10;
        x0 = ones(n,1);
        data = 'arglinb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 90
        n=10;
        x0 = ones(n,1);
        data = 'arglinc';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 91
        n=1000;
        x0 = ones(n,1);
        data = 'arwhead';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 92         %-1000，1001
        n=3;
        x0 = ones(n,1);
        data = 'bard';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 93
        n=1000;
        x0 = ones(n,1);
        data = 'bdexp';
        u = -10*ones(n,1);
        v = (10+1)*ones(n,1);
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 94
        n=1000;
        x0 = ones(n,1);
        data = 'bdqrtic';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 95
        n=6;
        x0 = [1 2 1 1 4 3]';
        data = 'biggs6';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 96
        n=3;
        x0 = [0 10 1]';
        data = 'box3';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 97
        n=2;
        x0 = [2 2]';
        data = 'brkmcc';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 98
        n=10;
        x0 = 1/2*ones(n,1);
        data = 'brownal';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 99
        n = 2;
        x0 = ones(n,1);
        data = 'brownbs';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 100
        n = 4;
        x0 = [25 5 -5 -1]';
        data = 'brownden';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 101
        n = 1000;
        x0 = ones(n,1);
        data = 'broydn7d';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 102
        ns = 499;
        n = 2*ns + 2;
        x0 = -2*ones(n,1);
        x0(1) = -3;
        x0(2) = -1;
        x0(3) = -3;
        x0(4) = -1;
        data = 'chainwoo';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 103
        n = 50;
        x0 = -ones(n,1);
        data = 'chnrosnb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 104
        n = 1000;
        x0 = ones(n,1);
        data = 'cosine';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 105
        n = 1000;
        x0 = 2*ones(n,1);
        x0(1) = 1;
        data = 'cragglvy';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 106
        n = 2;
        x0 = [-1.2 1]';
        data = 'cube';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 107
        n = 1000;
        xi = 1:n;
        x0 = (0.0001*xi/(n+1))';
        data = 'curly10';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 108
        n = 1000;
        x0 = 3*ones(n,1);
        data = 'dqdrtic';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 109
        n = 1000;
        x0 = 2*ones(n,1);
        data = 'dqrtic';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 110
        n = 1000;
        x0 = zeros(n,1);
        data = 'eg2';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 111
        n = 1000;
        x0 = 2*ones(n,1);
        data = 'engval1';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 112
        n = 3;
        x0 = [1 2 0]';
        data = 'engval2';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 113
        n = 50;
        x0 = -ones(n,1);
        data = 'errinros';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 114
        n = 10;
        x0 = ones(n,1);
        data = 'extrosnb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 115
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var/(h);
        data = 'fletcbv3';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 116
        n = 5;
        x0 = 506.2*ones(n,1);
        x0(1) = -506;
        data = 'genhumps';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 117
        n = 3;
        x0 = [5 2.5 0.15]';
        data = 'gulf';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 118
        n = 50;
        x0 = -3*ones(n,1);
        data = 'hilbertb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 119
        n = 2;
        x0 = [0 2]';
        data = 'himmelbh';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 120    %报错
        n = 4;
        x0 = [0.25 0.39 0.415 0.39]';
        data = 'kowosb';
        u = -10*ones(n,1);
        v = (10+1)*ones(n,1);

    case 121
        n = 3;
        x0 = [0.02 4000 250]';
        data = 'meyer3';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 122
        n = 1000;
        var = (1:n)';
        x0 = var;
        data = 'noncvxu2';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 123
        n = 5;
        x0 = [0.5 1.5 -1 0.01 0.02]';
        data = 'osbornea';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 124
        n = 100;
        x0 = 1/2*ones(n,1);
        data = 'penalty2';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 125
        n = 1000;
        x0 = 2*ones(n,1);
        data = 'quartc';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 126
        n = 2;
        x0 = [-1.2 1.0]';
        data = 'rosenbr';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 127
        n = 1000;
        var = (1:n)';
        scal = 12.0;
        scale = exp((var-1)*scal/(n-1));
        x0 = 1.0./scale;
        data = 'scosine';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 128
        n = 100;
        sc = 12;
        var1 = (0:n-1)';
        var2 = (1:n)';
        scale = exp(var1*sc/(n-1));
        x0 = 0.0001*(var2.*scale)/(n+1);
        data = 'scurly10';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 129
        n = 2;
        x0 = [4.712389 -1.0]';
        data = 'sineval';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 130
        n = 1000;
        x0 = 0.1*ones(n,1);
        data = 'sinquad';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 131
        n = 2;
        x0 = [1.0 0.1]';
        data = 'sisser';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 132
        n = 1000;
        x0 = ones(n,1);
        x0(1:2:n) = -1.2*x0(1:2:n);
        data = 'srosenbr';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 133
        n = 1000;
        x0 = ones(n,1);
        data = 'tridia';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 134
        n = 100;
        var = (1:n)';
        x0 = 1-(var/n);
        data = 'vardim';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 135
        n = 31;
        x0 = zeros(31,1);
        data = 'watson';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 136
        n = 1000;
        x0 = -ones(n,1);
        data = 'woods';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 137
        n = 2;
        x0 = [3 8]';
        data = 'zangwil2';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 138
        n = 1000;
        var = (1:n)';
        h = 1/(n+1);
        x0 = var*h;
        data = 'fletchbv';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 139
        n = 6;
        x0 = [0 0 1 1 1 1]';
        data = 'heart6ls';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 140
        n = 8;
        x0 = ones(n,1);
        data = 'heart8ls';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 141
        n = 10;
        x0 = ones(n,1);
        x0(1) = -4;
        x0(2) = -2;
        data = 'hilberta';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 142
        n = 4;
        x0 = [2.7 90 1500 10]';
        data = 'himmelbf';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 143
        n = 2;
        x0 = [-506.0 -506.2]';
        data = 'humps';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 144
        n = 31;
        x0 = [107.47 0.09203 0.908 102.4 0.1819 ...
            0.8181 97.44 0.284 0.716 96.3 ...
            0.3051 0.6949 93.99 0.3566 0.6434 ...
            89.72 0.468 0.532 83.71 0.6579 ...
            0.3421 78.31 0.8763 0.1237 886.37 ...
            910.01 922.52 926.46 935.56 952.83 ...
            975.73]';
        data = 'methanb8';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 145
        n = 2;
        x0 = [1.0e-30 1.0]';
        data = 'nasty';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 146
        n = 11;
        x0 = [1.3 0.65 0.65 0.7 0.6 ...
            3 5 7 2 4.5 5.5]';
        data = 'osborneb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);


    case 147
        n = 3;
        x0 = [0.6 -0.6 20]';
        data = 'yfitu';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 148
        n = 2;
        x0 = [1 1]';
        data = 'denschna';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 149
        n = 2;
        x0 = [1 1]';
        data = 'denschnb';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);

    case 150
        n = 3000;
        x0 = 2*ones(n,1);
        data = 'dixmaank';
        u = -1000*ones(n,1);
        v = (1000+1)*ones(n,1);


end
% data = 'Lev';	% select test function from Jones et al. test set
       
path(path,'jones');	% add path to directory with function 

% define bounds on variables (+-inf allowed)
%
% u: column vector of lower bounds
% v: column vector of upper bounds
% u(k)<v(k) is required
%

%%%%%%%%%%%%%%%2023/4/25修改
% [u,v,fglob] = defaults(data); 	% returns bounds used by Jones et al.
				% and known global optimum

% u = -1e+18*ones(n,1);
% v = (1e+18+1)*ones(n,1);
% u = -1000*ones(n,1);
% v = 10001*ones(n,1);
% u = -n*ones(n,1);
% v = (n+1)*ones(n,1);
% u = -1000*ones(n,1);
% v = (1000+1)*ones(n,1);


dimension=length(u)		% show dimension
%%%%%%%%%%%%%%%2023/4/25修改
% known_global_opt_value=fglob;	% show known global minimum value


% u=[-5,0]';v=[10,15]';    		% bra default bounds
% u=[-3,-2]';v=[3,2]';     		% cam default bounds
% u=[-2,-2]';v=[2,2]';     		% gpr default bounds
% u=[-10,-10]';v=[10,10]'; 		% shu default bounds
% u=[0,0,0]';v=[1,1,1]';   		% hm3 default bounds
% u=[0,0,0,0]';v=[10,10,10,10]';   	% sh5,sh7,s10 default bounds
% u=[0,0,0,0,0,0]';v=[1,1,1,1,1,1]';   	% hm6 default bounds
% modify the problem to be unconstrained by activating the next line
% u = -Inf*ones(size(u));v=-u; 


use_defaults=1;
% *** If you just want to use the default settings,
% *** you don't need to edit the rest of the file,
% *** If you are not satisfied with the results
% *** (this may happen especially when your box bounds are very wide),
% *** or if you want to play with the tuning parameters,
% *** set use_defaults=0, and modify the rest of the file 
% *** according to your curiosity or ingenuity

% if use_defaults, 
%   % easy to use version - all parameters preset
%   % defaults are being used, it suffices to call  
%   %%%%%%%%%%%%%%%%%% simple MCS call %%%%%%%%%%%%%%%%%%
%   [xbest,fbest,xmin,fmi,ncall,ncloc]=mcs(fcn,data,u,v);
%   % or, with less output,
%   % [xbest,fbest]=mcs(fcn,data,u,v);
% 
%   xbest	  		% best point found
%   fbest     		% best function value
%   fglob			% global minimum (known for test functions)
%   ncall	  		% number of function values used
%   ncloc	  		% number of function values in local searches
% 
%   % xmin	  	% columns are points in 'shopping basket'1
% 			% may be good alternative local minima
%   % fmi	  		% function values in 'shopping basket'
%   nbasket = length(fmi) % number of points in 'shopping basket'
% 
%   return;
% end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% change MCS settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flexible use version - all parameters may be modified

% 
% define amount of output printed
prt = 1;	% print level 
		% prt = 0: no output
		% prt = 1: # sweep, minimal nonempty level, # f-calls,
		%          best point and function value
		% prt > 1: in addition levels and function values of
		%          boxes containing the known global minimizers
		%          of a test function

% define global strategy 
%
% smax governs the relative amount of global versus local search. By
% increasing smax, more weight is given to global search.
% 
% Increasing or decreasing stop(1) increases or decreases the amount 
% of work allowed before concluding that nothing more is gained; 
% the default choice is quite conservative, and may try to locate
% a better point long after the global minimum has been found.
% stop(1)=5 works for many easier problems, too, with much fewer
% wasted function values.
% 
% Increasing nf allows to spend more time on boxes that have a chance 
% on their level to contain better points. This may be important for
% hard problems, or for problems with very wide bounds.
% 
% But in each case, it is unclear in advance what change would be most 
% beneficial to a particular problem. 
% We had very mixed experience; if you have many similar problems to 
% solve, the best thing to do is to experiment with a few problems to 
% find the best values, and to use these on the others. 
%
n = length(u);		% problem dimension
if (n > 500)
    smax = 5*n+10;		% number of levels used
    nf = n^2; 		% limit on number of f-calls 50*
else
    smax = 5*n+10;		% number of levels used
    nf = 50*n^2; 		% limit on number of f-calls 
end
stop(1) = 3*n;		% = m, integer defining stopping test 
stop(2) = -inf;		% = freach, function value to reach
			% if m>0, run until m sweeps without progress
			% if m=0, run until fbest<=freach
			% (or about nf function calls were used)

if 0, 	% known global optimum, for tests only
	% then the entries of stop have a different meaning
  stop(1) = 1.e-4;	% run until this relative error is achieved
  stop(2) = fglob;	% known global optimum value
  stop(3) = 1.e-10;	% stopping tolerance for tiny fglob
end;

% define initialization strategy
%
% for wide boxes, it is advisable (and for unbounded search regions
% strongly advisable) to define a customized initialization list
% that contains for each coordinate at least three reasonable values.
% Without such an initialization list, mcs makes default assumptions
% that roughly amount to estimating reasonable magnitudes from the 
% bounds and in case iinit=1 from assuming order unity if this is 
% within the bounds. 
%
% for a self-defined initialization list, the user should
% write an m-script file init0.m containing a matrix x0 with n
% rows and at least 3 columns and two n-vectors l and L 
% the ith column of x0 contains the initialization list
% values for the ith coordinate, their number is L(i), and
% x0(i,l(i)) is the ith coordinate of the initial point

iinit = 0;	% 0: simple initialization list
		%    (default for finite bounds;
		%     do not use this for very wide bounds)
		% 1: safeguarded initialization list
		%    (default for unbounded search regions)
		% 2: (5*u+v)/6, (u+v)/2, (u+5*v)/6
		% 3: initialization list with line searches
		% else: self-defined initialization list 
		%       stored in file init0.m

% parameters for local search
%
% A tiny gamma (default) gives a quite accurate but in higher 
% dimensions slow local search. Increasing gamma leads to less work 
% per local search but a less accurately localized minimizer
% 
% If it is known that the Hessian is sparse, providing the sparsity 
% pattern saves many function evaluations since the corresponding
% entries in the Hessian need not be estimated. The default pattern
% is a full matrix.
% 
local = 50;		% local = 0: no local search
			% else: maximal number of steps in local search
gamma = eps;		% acceptable relative accuracy for local search
hess = ones(n,n);	% sparsity pattern of Hessian



% defaults are not being used, use the full calling sequence
% (including at least the modified arguments)
%%%%%%%%%%%%%%%%%%%%%%% full MCS call %%%%%%%%%%%%%%%%%%%%%%
tic;
[xbest,fbest,xmin,fmi,ncall,ncloc]=...
  mcs(fcn,data,u,v,prt,smax,nf,stop,iinit,local,gamma,hess);
mcs_time = toc;
xbest	  		% best point found
fbest     		% best function value
ncall	  		% number of function values used
ncloc	  		% number of function values in local searches

% xmin	  		% columns are points in 'shopping basket'
			% may be good alternative local minima
% fmi	  		% function values in 'shopping basket'

solveNCPOther(test-68,2) = mcs_time;
solveNCPOther(test-68,1) = fbest;

% end


% nbasket = length(fmi) 	% number of points in 'shopping basket'
% MCS_result = 'The global minimum value found by MCS is %12.8f and the computational time of MCS is %6.4f s\n';
% fprintf(MCS_result,fbest,mcs_time);

