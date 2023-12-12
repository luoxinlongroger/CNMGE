

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% driverVRBBO.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file illustrates the use of VRBBO
% by showing how to minimize the function f(x):=||Ax-b||_p^e
% for p=2,e=1 in dimension n=300 from a random starting point.
% The goal assumed in the example is the reduction of the initial 
% objective function value by a factor of 0.01 with at most nfmax=200*n 
% function evaluations, to return when one of the two limiting 
% conditions is first satisfied.

clear;


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create problem definition

% define problem parameters (to be adapted to your problem)
n=input('>> choose dimension n = ');         % dimension
p=2;          % norm in objective function
e=1;          % exponent in objective function 

% create random matrix and right hand side for objective function
% (specific to the model problem; replace by whatever data you
% need to provide to your objectiv function)
A=rand(n)-0.5; 
b=-sum(A,2);

% create objective function f(x)=||Ax-b||_p^e
fun=@(x) norm(A*x-b,p).^e; 
% To solve your own problem, simply replace in the previous line 
% the expression after @(x) by your expression or function call. 
% Parameters in this expression are passed by value, hence are 
% fixed during minimization.

% start and stop info

x      = 2*rand(n,1)-1; % starting point
nfmax  = 200*n;         % stop after nfmax function evaluations
secmax = inf;           % stop after secmax seconds
finit  = fun(x);        % initial function value
accf   = 1e-3;          % stops stop when a point with 
fbest  = 0.01*finit;    % qf:= (f-fbest)/(finit-fbest) <= accf is found 

% For all CUTEst test problems, there have been a lot of efforts to
% get fbest by runing some competitive gradient-based solver and 
% gradient-free local and global solvers; for more details see the paper
% Here we use fbest = 0.01*finit for the above test problem 


% customization
% print level %  -1: nothing, 0: litte, >=1: more and more
prt=input(['>> print level (-1: nothing, 0: litte, >=1: more',...
           ' and more) = ']);
                
Tuning=0;   % 0: no tuning of VRBBO? (recommended at first)
            % 1: full tuning? (for specialists only)
            
fullinfo=1; % 0: pass stop and print criteria inside VRBBO
            % 1: pass stop and print criteria

% problem definition complete
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve the problem with VRBBO


if fullinfo
    % pass stop and print criteria
    % (indefinite run if no stopping criterion is given)
    st.secmax = secmax;       % (default: inf)
    st.nfmax  = nfmax;        % (default: inf)
    st.fbest  = fbest;        % (default: 0)
    st.accf   = accf;         % (default: 1e-4)
    st.finit  = finit;       
    st.prt    = prt;          % (default: -1)
else
   st = []; 
end

 

if Tuning
   % set tuning parameters (for specialists only)
   % given are defaults; these wouldn't need to be set 
  
    % maximum number of best points kept
    tune.msmax = 5;
    % maximum number of memory for L-BFGS
    tune.mqmax = 5; 
    % maximal number of multi-line searches in setScales
    tune.T0 = 2*tune.msmax;
    % maximal number of random trials
    tune.E= 10; 
    % scale for subspace direction
    tune.scSub = 0; 
    % scale for cumulative direction
    tune.scCum = 0; 
    % scale for heuristic direction
    tune.scCor = 0; 
    % select cumulative step type
    tune.cum = 1; 
    % upper bound for estimating cumulative step-size
    if tune.cum==2, tune.a = 3; end
    % minimium threshold for good impovement
    tune.Deltamin = 0;
    % initial threshold for good impovement
    tune.Deltamax = 1e-3; 
    % factor for finding delta
    tune.gammadelta = 1e6; 
    % factor for adjusting Deltamax
    tune.gammamax = 1e-3; 
    % factor for extrapolation test
    tune.gammaE = 4; 
    % factor for finding initial lam
    tune.gammalambda = 1e-6;
    % factor for angle condition
    tune.gammaw = eps;
    % factor for adjusting Delta
    tune.Q = 2; 
    % factor for angle condition
    tune.gammaangle = 1e-20; 
    % minimum norm of trial step
    tune.deltamin = 1e-4*sqrt(n);
    % maximum norm of trial step
    tune.deltamax = 0.1*sqrt(n);
    % minimium threshold for extrapolation step sizes
    tune.alphamin = 1e-50;
    % parameter for extrapolation test
    tune.gammamin = 1e-6;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify kind of algorithm by variable "alg":          %
    % 0: only random direction; R<n                          %
    % 1: only random direction; R >= n                       %
    % 2: no random subspace direction                        %
    % 3: random, random subspace, and cumulative directions  %
    % 4: no LBFGS direction                                  %
    % 5: all directions                                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tune.alg = 5;
else
    tune =[];
end
 
% call VRBBO 
[x,f,info] = VRBBO(fun,x,st,tune);
      
      
% display output

if prt>=0,
    disp('VRBBO completed silently'); 
    disp(' ');
    disp('display output');
    disp(' ');
    info               % progress report by VRBBO
    if prt>=1
       nstep=info.nstep
    end
    nfused=info.nf                  % number of function evaluations used
    secused=info.sec                % time used up to now
    if n<=100
       x'                           % best point found     
    end
    f                               % function value at xbest    
    qf=info.qf                      % target (for comparison)
    status = info.status
end;