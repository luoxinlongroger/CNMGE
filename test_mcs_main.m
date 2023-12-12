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
test = input('\n Please input the number of test problem:');
switch test
    case 1
        % molecular energy problem
        data = 'MoE';
        
    case 2
        % ackley_function
        data = 'Ack';
        
    case 3
        % levy_function
        data = 'Lev';
        
    case 4
        % schwefel_function
        data = 'Sch';
        
    case 5
        % rastrigin_function
        data = 'Ras';
        
    case 6
        % styblinski_tang_function
        data = 'StT';
        
    case 7
        % trid function
        data = 'Tri';
        
    case 8
        % sum_squares_function
        data = 'SuS';
        
    case 9
        % sphere_function
        data = 'Spe';
        
    case 10
        % rotated_hyper_ellipsoid_function
        data = 'RHE';
        
    case 11
        % zakharov_function
        data = 'Zak';
        
    case 12
        % dixon_price_function
        data = 'DiP';
        
    case 13
        % rosenbrock_function
        data = 'Ros';
        
    case 14
        % powell_function
        data = 'Pow';
        
    case 15
        % quartic_with_noise_function
        data = 'QWN';
        
    case 16
        % schubert_function
        data = 'Scu';
        
    case 17
        % raydan 1 function
        data = 'Ra1';
        
    case 18
        % raydan 2 function
        data = 'Ra2';
        
    case 19
        % extended tridiagonal-1 function
        data = 'ET1';
        
    case 20
        % extended quadratic penalty qp 1 function
        data = 'QP1';
        
    case 21
        % extended quadratic penalty qp 2 function
        data = 'QP2';
        
    case 22
        % quadratic qf2 function
        data = 'QF2';
        
    case 23
        % extended psc1 function
        data = 'PSC';
        
    case 24
        % extended bd1 function
        data = 'BD1';
        
    case 25
        % extended cliff function
        data = 'ECl';
        
    case 26
        % perturbed quadratic diagonal function
        data = 'Pqd';
        
    case 27
        % extended hiebert function -5 5
        data = 'Ehi';
        
    case 28
        % extended tet function
        data = 'TET';
        
    case 29
        % diagonal 1 function
        data = 'Di1';
        
    case 30
        % diagonal 3 function
        data = 'Di3';
        
    case 31
        % diagonal 5 function
        data = 'Di5';
        
    case 32
        % extended maratos function
        data = 'EMa';
        
    case 33
        % eg2 function
        data = 'EG2';
        
    case 34
        % sinquad function
        data = 'SIN';
        
    case 35
        % griewank_function
        data = 'Gri';
        
    case 36
        % levy13_function
        data = 'L13';
        
    case 37
        % hosaki_function
        data = 'Hos';
        
    case 38
        % beale_function
        data = 'Bea';
        
    case 39
        % easom_function
        data = 'Eas';
        
    case 40
        % price function
        data = 'Pri';
        
    case 41
        % branin_function
        data = 'Bra';
        
    case 42
        % trecanni_function
        data = 'Tre';
        
    case 43
        % booth_function
        data = 'Boo';
        
    case 44
        % matyas_function
        data = 'Mat';
        
    case 45
        % mccormick_function
        data = 'McC';
        
    case 46
        % power_sum_function
        data = 'PoS';
        
    case 47
        % colville_function
        data = 'Col';
        
    case 48
        % schaffer_function n.2
        data = 'SN2';
        
    case 49
        % bohachevsky_function
        data = 'Boh';
        
    case 50
        % three_hump_camel_function
        data = 'THC';
        
    case 51
        % six_hump_camel_function
        data = 'SHC';
        
    case 52
        % drop_wave_function
        data = 'DrW';
        
    case 53
        % perm_function
        data = 'Per';
        
    case 54
        % hartmann_3_d_dimensional_function
        data = 'Hm3';
        
    case 55
        % trefethen_4_function
        data = 'Tr4';
        
    case 56
        % zettl_function
        data = 'Zet';
        
    case 57
        % exp2_function
        data = 'Ex2';
        
    case 58
        % hansen_function
        data = 'Han';
        
    case 59
        % schaffer_function_n_4
        data = 'SN4';
        
    case 60
        % holder_table_function
        data = 'HoT';
        
    case 61
        % gramacy_lee_function
        data = 'GrL';
        
    case 62
        % eggholder_function
        data = 'Egg';
        
    case 63
        % michalewicz_function
        data = 'Mic';
        
    case 64
        % box_betts_exponential_quadratic_sum_function
        data = 'BBE';
        
    case 65
        % cross_in_tray_function
        data = 'CiT';
        
    case 66
        % himmelblau_function
        data = 'Him';
        
    case 67
        % forrester_function
        data = 'For';
        
    case 68
        % goldstein-price function
        data = 'GPr';
end
% data = 'Lev';	% select test function from Jones et al. test set
       
path(path,'jones');	% add path to directory with function 

% define bounds on variables (+-inf allowed)
%
% u: column vector of lower bounds
% v: column vector of upper bounds
% u(k)<v(k) is required
%
[u,v,fglob] = defaults(data); 	% returns bounds used by Jones et al.
				% and known global optimum
dimension=length(u)		% show dimension
known_global_opt_value=fglob;	% show known global minimum value
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
nbasket = length(fmi) 	% number of points in 'shopping basket'
MCS_result = 'The global minimum value found by MCS is %12.8f and the computational time of MCS is %6.4f s\n';
fprintf(MCS_result,fbest,mcs_time);