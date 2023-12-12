% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% This code is revised by Xin-long Luo, on Augus 14-22, 2022, in order to
% process the case without the analytical gradient of the objective
% function and use automatic diffferentiation with the reverse mode to
% compute its gradient. 
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% senzhang@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% This solver detetermines the global minmum of unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 1
% and quasi-genetic ecolution. This program runs all programs successively 
% in desired way.
%
% This is an entrance subroutine when we use this software to obtain the
% global optimal solution of unconstrained optimization problem.
%
% Have a look at "Xin-long Luo, Hang Xiao, Sen Zhang 
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
% func: it denotes the test problems.
%
% lb: it denotes the lower boundary of the variable.
%
% ub: it denotes the upper boundary of the variable.
%
% Output: 
% There is no output parameter. But we output the global optimal solution 
% of the problem, the number of local optimal solutions obtained by CNMGE
% and the computational time.
%
% We call gcp('nocreate') to get Get current parallel pool.
% poolobj = gcp('nocreate');% If no pool,  create new one.
% if isempty(poolobj)
%     poolsize = 0;
%     CoreNum=4; 
%     parpool(CoreNum);
% else
%     poolsize = poolobj.NumWorkers;
%     disp('Already initialized');
% end

clear

test = input('\n Please input the number of test problem:');
switch test
    case 1 
        % molecular energy problem 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function 
        x = optimvar('x',n,1);
        k = 4.141720682;
        r = 10.60099896;
        p = 1 + cos(3*x);         

        % Generate an n-dimensional vector [-1,+1,-1,+1,...,-1,+1].
        allones = ones(n,1);
        odds = 1:2:n;
        oddeven = allones;
        oddeven(odds) = allones(odds) - 2*allones(odds);
        q = oddeven./sqrt(r - k*cos(x));        
        objfun = sum(p) + sum(q);  
        
        lb = -inf;
        ub = inf;
        
    case 2  
        % ackley_function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function 
        x = optimvar('x',n,1);
        d = 1000;
        c = 2*pi;
        b = 0.2;
        a = 20;
        sum1 = x'*x;
        cosy = cos(c*x);
        sum2 = sum(cosy);
        term1 = -a * exp(-b*sqrt(sum1/d));
        term2 = -exp(sum2/d);
        objfun = term1 + term2 + a + exp(1);
        
        lb = -inf;
        ub = inf;
        
    case 3
        % levy_function
        n = 1000;
        x0 = 2*ones(n,1);
        % Define the objective function 
        x = optimvar('x',n,1);

        w = 1 + (x-1)/4;
        term1 = (sin(pi*w(1)))^2;
        term3 = (w(n)-1)^2 * (1+(sin(2*pi*w(n)))^2);

        wns1 = w(1:(n-1)); 
        wns1s1 = wns1 - 1;
        wns1s1sq = wns1s1.^(2);
        vectwns1 = sin(pi*wns1+1);
        vectwns1sq = vectwns1.^(2);
        term2vec = wns1s1sq.*(1+10*vectwns1sq);
        term2 = sum(term2vec);

        objfun = term1 + term2 + term3;

        lb = -inf;
        ub = inf;
        
    case 4  
        % schwefel_function
        n = 1000;
        x0 = 200*ones(n,1);
        % Define the objective function 
        x = optimvar('x',n,1);
        xsqr = x.^2;
        absx = sqrt(xsqr);      
        sqrtabsx = sqrt(absx);
        sinvec = sin(sqrtabsx);
        vectmid = x.*sinvec;

        objfun = 418.9829*n - sum(vectmid);

        lb = -500;
        ub = 500;

    case 5
        % rastrigin_function 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        xsq = x.^2;
        cosvec = cos(2*pi*x);
        vectm = xsq - 10*cosvec;

        objfun = 10*n + sum(vectm);

        lb = -inf;
        ub = inf;
        
    case 6
        % styblinski_tang_function 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        vectvar = x.^(4) - 16*x.^(2) + 5*x;
        objfun = sum(0.5*vectvar);

        lb = -inf;
        ub = inf;
        
    case 7
         % trid function
         n = 1000;
         x0 = ones(n,1);
         % Define the objective function
         x = optimvar('x',n,1);
         x1tns1 = x(1:(n-1));
         x2tn = x(2:n);  
         vectvar1 = x - 1;            
         term1 = dot(vectvar1,vectvar1);        
         term2 = dot(x2tn,x1tns1);          
         objfun =  term1 - term2;               

         lb = -inf;
         ub = inf;
         
    case 8
        % sum_squares_function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        xsqr = x.^2;
        vectvar = indexvart.*xsqr;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
     case 9
        % sphere_function 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        objfun = dot(x,x);

        lb = -inf;
        ub = inf;
        
     case 10
        % rotated_hyper_ellipsoid_function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        xsqr = x.^2;
        invindex = n:(-1):1;
        invindext = invindex';
        vectvar = invindext.*xsqr;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
     case 11
        % zakharov_function 
        n = 1000;
        x0 = (5e-5)*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n; 
        indexvart = indexvar';
        xsqr = x.^2;        
        sumfun1 = sum(xsqr);
        vectvar = 0.5*(indexvart.*x);
        sumfun2 = sum(vectvar);
        objfun = sumfun1 + sumfun2^2 + sumfun2^4;

        lb = -inf;
        ub = inf;
        
     case 12
        % dixon_price_function 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        term1 = (x(1)-1)^2;
        x2tn = x(2:n);
        x1tns1 = x(1:(n-1));        
        index2tn = 2:n;
        index2tnt = index2tn';
        vectvar1 = 2*(x2tn.^2) - x1tns1;
        vectvar1sqr = vectvar1.^2;
        vectvar2 = index2tnt.*vectvar1sqr;
        objfun = term1 + sum(vectvar2);

        lb = -inf;
        ub = inf;
        
    case 13
        % rosenbrock_function
        n = 1000;
        x0 = 2*ones(n,1);
        % Define the objective function 
        x = optimvar('x',n,1);
        x1tns1 = x(1:n-1);
        x2tn = x(2:n);
        vectvar1 = x2tn - x1tns1.^2;
        vectvar1sqr = vectvar1.^2;
        vectvar2 = x1tns1 - 1;
        vectvar2sqr = vectvar2.^2;
        objfun = 100*sum(vectvar1sqr) + sum(vectvar2sqr);

        lb = -inf;
        ub = inf;
        
    case 14
        % powell_function
        n = 1000;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        index0mod4 = 4:4:n;
        index1mod4 = 1:4:n;
        index2mod4 = 2:4:n;
        index3mod4 = 3:4:n;
        vectvar1 = x(index1mod4) + 10*x(index2mod4);
        vectvar1sqr = vectvar1.^2;

        vectvar2 = x(index3mod4) - x(index0mod4);
        vectvar2sqr = vectvar2.^2;

        vectvar3 = x(index2mod4) - 2*x(index3mod4);
        vectvar3p4 = vectvar3.^4;

        vectvar4 = x(index1mod4) - x(index0mod4);
        vectvar4p4 = vectvar4.^4;

        objfun = sum(vectvar1sqr) + 5*sum(vectvar2sqr)...
            + sum(vectvar3p4) + 10*sum(vectvar4p4);

        lb = -inf;
        ub = inf;
        
    case 15
         % quartic_with_noise_function 
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        noise = 0.5;
        xquart = x.^4;
        objfun = sum(xquart) + noise;

        lb = -inf;
        ub = inf;
        
    case 16
        % schubert_function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        sin2xp1 = sin(2*x+1);
        sin3xp2 = sin(3*x+2);
        sin4xp3 = sin(4*x+3);
        sin5xp4 = sin(5*x+4);
        sin6xp5 = sin(6*x+5);
        sumvec = sin2xp1 + 2*sin3xp2 + 3*sin4xp3 + 4*sin5xp4 + 5*sin6xp5;
        objfun = -sum(sumvec);

        lb = -10;
        ub = 10;
        
    case 17
        % raydan 1 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        vectvar = indexvart.*(exp(x) - x);
        objfun = sum(vectvar)/10;

        lb = -inf;
        ub = inf;
        
    case 18
        % raydan 2 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        vectvar = exp(x) - x;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 19
        % extended tridiagonal-1 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        vectvar1 = x(oddindex) + x(evenindex) - 1;
        vectvar1sqr = vectvar1.^2;
        vectvar2 = x(oddindex) - x(evenindex) + 1;
        vectvar2qua = vectvar2.^4;
        objfun = sum(vectvar1sqr + vectvar2qua);

        lb = -inf;
        ub = inf;
        
     case 20
        % extended quadratic penalty qp 1 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        xsqr = x.^2;
        sumx2 = sum(xsqr) - 0.5;
        x1tn = x(1:(n-1));
        xsqr2 = x1tn.^2;
        xsqrs2 = xsqr2 - 2;
        xquad = xsqrs2.^2;
        objfun = sum(xquad) + sumx2^2;

        lb = -inf;
        ub = inf;
        
    case 21
        % extended quadratic penalty qp 2 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        xsqr = x.^2;
        sumx2 = sum(xsqr) - 100;

        x1tn = x(1:(n-1));
        xsqr2 = x1tn.^2;
        vectvar = xsqr2 - sin(x1tn);
        vectvarsqr = vectvar.^2;
        objfun = sum(vectvarsqr) + sumx2^2;

        lb = -inf;
        ub = inf;
        
    case 22
        % quadratic qf2 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        xsqr = x.^2;
        vectvar = xsqr - 1;
        vectvarsqr = vectvar.^2;
        vectvarsqrmid = indexvart.*vectvarsqr;
        objfun = 0.5*sum(vectvarsqrmid) - x(n);

        lb = -inf;
        ub = inf;
        
    case 23
        % extended psc1 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;        
        xodd = x(oddindex);
        xeven = x(evenindex);
        xoddsqr = xodd.^2;
        xevensqr = xeven.^2;
        xoddmxeven = xodd.*xeven;
        vectsqr = xoddsqr + xevensqr + xoddmxeven;
        vectquad = vectsqr.^2;

        sinxodd = sin(xodd);
        sinxoddsqr = sinxodd.^2;
        cosxeven = cos(xeven);
        cosxevensqr = cosxeven.^2;

        vectvar = vectquad + sinxoddsqr + cosxevensqr;      
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 24
        % extended bd1 function
        n = 1000;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        xodd = x(oddindex);
        xeven = x(evenindex);
        xoddsqr = xodd.^2;
        xevensqr = xeven.^2;

        vectvar1 = xoddsqr + xevensqr - 2;
        vectvar1sqr = vectvar1.^2;

        vectvar2 = exp(xodd -1) - xeven;
        vectvar2sqr = vectvar2.^2;

        objfun = sum(vectvar1sqr + vectvar2sqr);

        lb = -inf;
        ub = inf;
        
    case 25
        % extended cliff function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        xodd = x(oddindex);
        xeven = x(evenindex);

        vectvar1 = (xodd - 3)/100;
        vectvar1sqr = vectvar1.^2;

        vectvar2 = xodd - xeven;
        vectvar3 = exp(20*vectvar2);

        vectvar = vectvar1sqr - vectvar2 + vectvar3;        
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 26
        % perturbed quadratic diagonal function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        xsqr = x.^2;
        vectvar = indexvart.*xsqr;
        sumx = sum(x);

        objfun = (sum(vectvar))/100 + sumx^2;

        lb = -inf;
        ub = inf;
         
    case 27
        % extended hiebert function
        n = 1000;
        x0 = 2000*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        xodd = x(oddindex);
        xeven = x(evenindex);

        vectvar1 = xodd - 10;
        vectvar1sqr = vectvar1.^2;

        vectvar2 = xodd.*xeven - 50000;
        vectvar2sqr = vectvar2.^2;

        vectvar = vectvar1sqr + vectvar2sqr;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 28
        % extended tet function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        xodd = x(oddindex);
        xeven = x(evenindex);

        vectvar1 = exp(xodd+3*xeven-0.1);
        vectvar2 = exp(xodd-3*xeven-0.1);
        vectvar3 = exp(-xodd-0.1);

        vectvar = vectvar1 + vectvar2 + vectvar3;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 29
        % diagonal 1 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        vectvar1 = exp(x);
        vectvar2 = indexvart.*x;
        vectvar = vectvar1 - vectvar2;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
   
    case 30
        % diagonal 3 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        vectvar1 = exp(x);
        vectvar2 = indexvart.*sin(x);
        vectvar = vectvar1 - vectvar2;
        objfun = sum(vectvar);

        lb = -inf;
        ub = inf;
        
    case 31
        % diagonal 5 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        vectvar = exp(x) + exp(-x);
        vectvarlog = log(vectvar);
        objfun = sum(vectvarlog);

        lb = -inf;
        ub = inf;
        
    case 32
        % extended maratos function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        oddindex = 1:2:n;
        evenindex = 2:2:n;
        xodd = x(oddindex);
        xeven = x(evenindex);

        xoddsqr = xodd.^2;
        xevensqr = xeven.^2;
        vectvar = xoddsqr + xevensqr - 1;
        vectvarsqr = vectvar.^2;

        objfun = sum(xodd+100*vectvarsqr);

        lb = -inf;
        ub = inf;
        
    case 33
        % eg2 function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        x1tns1 = x(1:(n-1));
        x1tns1sqr = x1tns1.^2;
        vectvar = x(1) + x1tns1sqr - 1;
        sinvectvar = sin(vectvar);
        objfun = sum(sinvectvar) + 0.5*sin(x(n)^2);

        lb = -inf;
        ub = inf;
        
    case 34
        % sinquad function
        n = 1000;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        vectvar1 = (x(1)-1)^4 + (x(n)^2-x(1)^2)^2;
        x2tns1 = x(2:(n-1));
        x2tns1sqr = x2tns1.^2;
        sinx2tns1 = sin(x2tns1-x(n));
        vectvar2 = sinx2tns1 - x(1)^2 + x2tns1sqr;
        vectvar2sqr = vectvar2.^2;
        objfun = vectvar1 + sum(vectvar2sqr);

        lb = -inf;
        ub = inf;
        
    case 35  
         % griewank_function
        n = 10;
        x0 = 1*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        indexvart = indexvar';
        indexvarsqrt = sqrt(indexvart);
        xsqr = x.^2;
        sumfun = sum(xsqr)/4000;
        vectvar = x./indexvarsqrt;
        cosvectvar = cos(vectvar);
        prodfun = prod(cosvectvar);
        objfun = sumfun - prodfun + 1;

        lb = -inf;
        ub = inf;
        
    case 36  
         % levy13_function
        n = 2;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        
        term1 = (sin(3*pi*x1))^2;
        term2 = ((x1-1)^2)*(1+(sin(3*pi*x2))^2);
        term3 = ((x2-1)^2)*(1+(sin(2*pi*x2))^2);
        
        objfun = term1 + term2 + term3;

        lb = -inf;
        ub = inf;

    case 37
        % hosaki_function
        n = 2;
        x0 = -2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        term1 = 1-8*x(1)+7*x(1)^2-(7*x(1)^3)/3+(x(1)^4)/4;
        term2 = x(2)^2;
        term3 = exp(-x(2));

        objfun = term1*term2*term3;
        
        lb = -inf;
        ub = inf;

    case 38
        % beale_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        
        term1 = (1.5 - x1 + x1*x2)^2;
        term2 = (2.25 - x1 + x1*x2^2)^2;
        term3 = (2.625 - x1 + x1*x2^3)^2;
        
        objfun = term1 + term2 + term3;

        lb = -inf;
        ub = inf;
        
    case 39
        % easom_function
        n = 2;
        x0 = -1*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);

        lb = -inf;
        ub = inf;
        
    case 40
        % price function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = (2*x(1)^3*x(2)-x(2)^3)^2+(6*x(1)-x(2)^2+x(2))^2;

        lb = -inf;
        ub = inf;
        
    case 41
        % branin_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        t = 1 / (8*pi);
        s = 10;
        r = 6;
        c = 5/pi;
        b = 5.1/(4*pi^2);
        a = 1;
        
        term1 = a*(x2 - b*x1^2 + c*x1 - r)^2;
        term2 = s*(1-t)*cos(x1);
        
        objfun = term1 + term2 + s;

        lb = -inf;
        ub = inf;
        
    case 42
        % trecanni_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = x(1)^4 + 4*x(1)^3 + 4*x(1)^2 + x(2)^2;

        lb = -inf;
        ub = inf;

    case 43
        % booth_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;

        lb = -inf;
        ub = inf;
        
    case 44
        % matyas_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = 0.26*(x(1)^2+x(2)^2)-0.48*x(1)*x(2);

        lb = -inf;
        ub = inf;
        
    case 45
        % mccormick_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = sin(x(1)+x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1;

        lb = -3;
        ub = 4;
        
    case 46
        % power_sum_function 
        n = 4;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        b = [8;18;44;114];

        outer = 0;        
        for ii = 1:n
            xpii = x.^ii;
            inner = sum(xpii);
            outer = outer + (inner-b(ii))^2;
        end
        objfun = outer;

        lb = -inf;
        ub = inf;
        
    case 47
        % colville_function
        n = 4;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        
        term1 = 100*((x1^2-x2)^2);
        term2 = (x1-1)^2;
        term3 = (x3-1)^2;
        term4 = 90*((x3^2-x4)^2);
        term5 = 10.1*((x2-1)^2 + (x4-1)^2);
        term6 = 19.8*(x2-1)*(x4-1);

        objfun = term1 + term2 + term3 + term4 + term5 + term6;

        lb = -inf;
        ub = inf;
        
    case 48
        % schaffer_function_n_2
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        numerator = sin(sqrt(x(1)^2+x(2)^2))^2-0.5;
        denominator = (1+0.001*(x(1)^2+x(2)^2))^2;
        objfun = 0.5 + numerator/denominator;

        lb = -inf;
        ub = inf;

    case 49
        % bohachevsky_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;

        lb = -inf;
        ub = inf;
        
    case 50
        % three_hump_camel_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = 2*x(1)^2-1.05*x(1)^4+x(1)^6/6+x(1)*x(2)+x(2)^2;

        lb = -inf;
        ub = inf;
        
    case 51
        % six_hump_camel_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        objfun = (4-2.1*x(1)^2+x(1)^4/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;

        lb = -inf;
        ub = inf;
        
    case 52
        % drop_wave_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = -(1+cos(12*sqrt(x(1)^2+x(2)^2)))/(0.5*(x(1)^2+x(2)^2)+2);

        lb = -inf;
        ub = inf;
        
    case 53
        % perm_function
        n = 4;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        b = 10;
        indexvar = 1:n;
        indexvart = indexvar';
        outer = 0;
        for ii = 1:n
            idpii = indexvart.^ii;
            idpiiaddb = idpii + b;
            xdid = x./indexvart;
            xdidpii = xdid.^ii;
            xdipiisub1 = xdidpii - 1;
            inner = sum(idpiiaddb.*xdipiisub1);
            outer = outer + inner^2;
        end

        objfun = outer;

        lb = -inf;
        ub = inf;
        
     case 54
        % hartmann_3_d_dimensional_function
        n = 3;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        alpha = [1.0;1.2;3.0;3.2];

        A = [3.0, 10, 30;
             0.1, 10, 35;
             3.0, 10, 30;
             0.1, 10, 35];

        B = [3689, 1170, 2673;
             4699, 4387, 7470;
             1091, 8732, 5547;
             381, 5743, 8828];
        P = 10^(-4)*B;

        outer = 0;
        jj = 1:3;
        for ii = 1:4
            vectvar1 = (x-P(ii,:)').^2;
            vectvar2 = A(ii,:)';
            inner = sum(vectvar1.*vectvar2);
            outer = outer + alpha(ii)*exp(-inner);
        end
        objfun = -outer;        

        lb = -inf;
        ub = inf;
        
     case 55
        % trefethen_4_function
        n = 2;
        x0 = 0.5*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        term1 = exp(sin(50*x(1)));
        term2 = sin(60*exp(x(2)));
        term3 = sin(70*sin(x(1)));
        term4 = sin(sin(80*x(2)));
        term5 = -sin(10*(x(1)+x(2)));
        term6 = 0.25*(x(1)^2+x(2)^2);

        objfun = term1 + term2 + term3 + term4 + term5 + term6;        

        lb = -inf;
        ub = inf;
        
     case 56
        % zettl_function
        n = 2;
        x0 = 0.1*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        objfun = (x(1)^2+x(2)^2-2*x(1))^2+0.25*x(1);

        lb = -inf;
        ub = inf;
        
     case 57
        % exp2_function
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 0:9;
        term1 = exp(-indexvar*(x(1)/10));
        term2 = -5*exp(-indexvar*(x(2)/10));
        term3 = -exp(-indexvar/10);
        term4 = 5*exp(-indexvar);
        vectvar = term1 + term2 + term3 + term4;
        vectvarsqr = vectvar.^2;
        objfun = sum(vectvarsqr);

        lb = -inf;
        ub = inf;
        
     case 58
        % hansen_function
        n = 2;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:5;
        cosindexvar1 = cos((indexvar+1)*x(2)+indexvar);
        vectvar1 = indexvar.*cosindexvar1;
        temp1 = sum(vectvar1);

        cosindexvar2 = cos((indexvar-1)*x(1)+indexvar);
        vectvar2 = indexvar.*cosindexvar2;
        temp2 = sum(vectvar2);
        objfun = temp1*temp2;

        lb = -inf;
        ub = inf;
        
     case 59
        % schaffer_function_n_4
        n = 2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        varx = x1^2 - x2^2;        
        absvarx = norm(varx);
        
        temp1 = (cos(sin(absvarx)))^2 - 0.5;
        temp2 = (1 + 0.001*dot(x,x))^2;
        
        objfun = 0.5 + temp1/temp2;

        lb = -inf;
        ub = inf;
        
     case 60
        % holder_table_function
        n = 2;
        x0 = 10*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        nx = norm(x,2);
        nxs1 = 1 - nx/pi;
        nxs1sqr = nxs1^2;
        absnxs1 = sqrt(nxs1sqr);
        
        temp1 = sin(x1)*cos(x2);
        temp2 = exp(absnxs1);

        pd12 = temp1*temp2;
        pd12sqr = pd12^2;
        abspd12 = sqrt(pd12sqr);
        
        objfun = -abspd12;

        lb = -inf;
        ub = inf;
        
     case 61
        % gramacy_lee_function
        n = 1;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = sin(10*pi*x(1))/(2*x(1))+(x(1)-1)^4;

        lb = 0.5;
        ub = 3;
        
     case 62
        % eggholder_function
        n = 2;
        x0 = 500*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);
        t1 = x2+x1/2+47;
        t1sqr = t1^2;
        abst1 = sqrt(t1sqr); 

        t2 = x1-(x2+47);
        t2sqr = t2^2;
        abst2 = sqrt(t2sqr); 
        
        term1 = -(x2+47)*sin(sqrt(abst1));
        term2 = -x1*sin(sqrt(abst2));
        
        objfun = term1 + term2;

        lb = -inf;
        ub = inf;

     case 63
        % michalewicz_function 
        n = 2;
        x0 = 2*ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        indexvar = 1:n;
        indexvart = indexvar';
        xsqr = x.^2;
        sinx = sin(x);
        idmxsqr = (indexvart.*xsqr)/pi;
        sinix = sin(idmxsqr);
        sinixp20 = sinix.^20;
        objfun = -sum(sinx.*sinixp20);

        lb = -inf;
        ub = inf;
        
    case 64
        % box_betts_exponential_quadratic_sum_function
        n = 3;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        indexvar = 1:n;
        term1 = exp(-0.1*indexvar*x(1));
        term2 = exp(-0.1*indexvar*x(2));
        term3 = (exp(-0.1*indexvar)-exp(-indexvar))*x(3);
        
        vectvar = term1 - term2 - term3;
        vectvarsqr = vectvar.^2;
        objfun = sum(vectvarsqr);

        lb = -inf;
        ub = inf;
        
    case 65
        % cross_in_tray_function
        n=2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);

        dx = norm(x,2);
        t1 = 100 - dx/pi;
        t1sqr = t1^2; 
        abst1 = sqrt(t1sqr);
        
        fact1 = sin(x1)*sin(x2);
        fact2 = exp(abst1);

        t2 = fact1*fact2;
        t2sqr = t2^2;
        abst2 = sqrt(t2sqr);
        
        objfun = -0.0001*(abst2+1)^0.1;

        lb = -10;
        ub = 10;
        
    case 66
        % himmelblau_function
        n=2;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        objfun = (x1^2 + x2 -11)^2 + (x1 + x2^2 -7)^2;

        lb = -inf;
        ub = inf;
        
    case 67
        % forrester_function
        n=1;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = (6*x(1)-2)^2*sin(12*x(1)-4);

        lb = -inf;
        ub = inf;
        
    case 68
        % goldstein-price function
        n = 2;
        x0 = [0;1];
        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);
        term1 = 1+(x1+x2+1)^2*(19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2);
        term2 = 30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2);
        objfun = term1*term2;   

        lb = -2;
        ub = 2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%cute
    case 69
        %allinitu       12868194
        n = 4;
        x0 = ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        objfun = 	x(3)-1 + x(1)^2+ x(2)^2 + (x(3)+x(4))^2 + sin(x(3))^2 + ...
        x(1)^2*x(2)^2 + x(4)-3 + ...
    	sin(x(3))^2 + (x(4)-1)^2 + (x(2)^2)^2 +... 
        (x(3)^2 + (x(4)+x(1))^2)^2 + (x(1)-4 + ...
        sin(x(4))^2 + x(2)^2*x(3)^2)^2 + sin(x(4))^4;

        lb = -inf;
        ub = inf;

%     case 70
%         %bdexp
%         n = 500;
%         x0 = ones(n,1);
%         % Define the objective function
%         x = optimvar('x',n,1);
%         
%         term1 = x(1:1:n-2);
%         term2 = x(2:1:n-1);
%         term3 = x(3:1:n);
% 
%         objfun = sum((term1+term2).*exp((term1+term2).*(-term3)));
%         lb = -inf;
%         ub = inf;

    case 71 
        %cliff  12868197
        n = 2;
        x0 = [0;-1];
        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);
        
        objfun = (0.01*x1-0.03)^2 -x1+x2+exp(20*(x1-x2));
        
        lb = -inf;
        ub = inf;

    case 72 
        %dixon3dq   12868200
        n = 10;
        x0 = -ones(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(2:1:n-1);
        term2 = x(3:1:n);
        
        objfun = (x(1)-1.0)^2 + ...
                sum((term1-term2).^2) + (x(10)-1.0)^2;
        
        lb = -inf;
        ub = inf;

    case 73
        %edensch        12868204
        n = 2000;
        x0 = zeros(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:1:n-1);
        term2 = x(2:1:n);
        
        objfun =  sum( (term1-2).^4 + (term1.*term2- ...
            2*term2).^2 + (term2+1).^2 ) + 16;
        
        lb = -inf;
        ub = inf;

    case 74
        %fletchcr   12868280 
        n = 100;
        x0 = zeros(n,1);
        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:1:n-1);
        term2 = x(2:1:n);
        
        objfun =  sum(100*(term2-term1+1-term1.^2).^2);
        
        lb = -inf;
        ub = inf;

    case 75
        %genrose    12868285 
        n = 500;
        x0 = ones(n,1)/501;
        
        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(2:1:n);
        term2 = x(1:1:n-1);
        
        objfun =  	1.0 +...
	                sum (100*(term1-term2.^2).^2) +...
	                sum ((term1-1.0).^2);
        
        lb = -inf;
        ub = inf;

    case 76
        %hairy  12868287
        n = 2;
        x0 = [-5;-7];
        hlength = 30;
        cslope = 100;
        
        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);

        objfun = sin(7*x1)^2*cos(7*x2)^2*hlength + ...
                 cslope*sqrt(0.01+(x1-x2)^2) + ...
                 cslope*sqrt(0.01+x1^2);
        
        lb = -inf;
        ub = inf;

    case 77
        %himmelbb   12868288
        n = 2;
        x0 = [-1.2;1];

        % Define the objective function
        x = optimvar('x',n,1);

        objfun = (x(1)*x(2)*(1-x(1))*(1-x(2)-x(1)*(1-x(1)^5)))^2;

        lb = -inf;
        ub = inf;

    case 78
        %himmelbg   12868326
        n = 2;
        x0 = [0.5;0.5];

        % Define the objective function
        x = optimvar('x',n,1);

        objfun = exp(-x(1)-x(2))*(2*x(1)^2+3*x(2)^2);

        lb = -inf;
        ub = inf;

    case 79
        %indef  12868331
        n = 1000;
        alpha = 0.5;
        x0 = 1/1001*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        xn = x(n);
        term = x(2:1:n-1);

        objfun = sum(x) +...
	            sum (alpha*cos(2*term-xn-x1));

        lb = -inf;
        ub = inf;

    case 80
        %jensmp 12868341
        n = 2;
        x0 = [0.3;0.4];

        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);
        indexvar =  1:10;
        indexvart = indexvar';
        vectvar1 = x1*indexvart;
        vectvar2 = x2*indexvart;

        objfun = sum((2+2*indexvart-(exp(vectvar1)+exp(vectvar2))).^2);


%         for i=1:10
%             objfun =  objfun + (2+2*i-(exp(i*x1)+exp(i*x2)))^2;
%         end

        lb = -inf;
        ub = inf;

    case 81
        %liarwhd    12868354
        n = 1000;
        x0 = 4*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);

        objfun = sum(4*((-x1+x.^2).^2)) + ...
                 sum((x-1.0).^2);

        lb = -inf;
        ub = inf;

    case 82
        %loghairy   12868358 
        n = 2;
        hlength = 30;
        cslope = 100;
        x0 = [-500;-700];

        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);

        objfun = log( (100+sin(7*x1)^2*cos(7*x2)^2*hlength +...
                cslope*sqrt(0.01+(x1-x2)^2) + ...
                cslope*sqrt(0.01+x1^2))/100);

        lb = -inf;
        ub = inf;  

    case 83
        %maratosb   12868361
        n = 2;
        invp = 0.000001;
        x0 = [0;0];

        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);

        objfun = x1+(x1^2+x2^2-1)^2/invp;

        lb = -inf;
        ub = inf;

    case 84
        %mexhat 12868363
        n = 2;
        x0 = [0.86;0.72];
        p = 10000;

        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);

        objfun = - 2*(x1-1)^2 + ...
                    p*(-0.02 + (x2-x1^2)^2/p + (x1-1)^2)^2;

        lb = -inf;
        ub = inf;  


    case 85
        %nondia 12957752    原题为1w维
        n = 1000;
        x0 = -1*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term = x(1:n-1);

        objfun =  	(x(1)-1)^2 +... 
                    sum(100*(x(1)-term.^2).^2);

        lb = -inf;
        ub = inf;

    case 86
        %nondquar   原题为1w维  12957833 
        n = 1000;
        x0 = ones(n,1);

        for i=1:n
            if(mod(i,2) == 0)
                x0 = -x0;
            end 
        end 
        

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-2);
        term2 = x(2:n-1);
        xN = x(n);

        objfun =  sum((term1+term2+xN).^4 + (x(1)-x(2))^2 + (x(n-1)+xN)^2);

        lb = -inf;
        ub = inf;      

    case 87
        %penalty1   12868384
        n = 1000;
        a = 10^-5;
        x0 = 1:n;
        x0 = x0';

        % Define the objective function
        x = optimvar('x',n,1);

        objfun =  sum(a*((x-1).^2)) + (sum (x.^2 - 1/4))^2;

        lb = -inf;
        ub = inf;  

    case 88
        %power  12868387
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        indexvar =  1:n;
        indexvart = indexvar';

        objfun =  sum((x.*indexvart).^2);

        lb = -inf;
        ub = inf; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 case 89
     %arglinb   %%%%%正确
     n = 10;
     m = 20;
     x0 = ones(n,1);

     % Define the objective function
     x = optimvar('x',n,1);

     indexvarM = 1:m;
     indexvarMt = repmat(indexvarM,n,1);
     indexvarN = 1:n;

     objfun = sum(((x'.*indexvarN)*indexvarMt-1).^2);

     lb = -inf;
     ub = inf;

    case 90
        %arglinc    %%%%%%
        n = 10;
        m = 20;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        indexvarN = 1:10;
        indexvarN(1) = 0;
        indexvarN(10) = 0;
        indexvarM = 1:m-1;
        indexvarM = indexvarM-1;
        indexvarMt = repmat(indexvarM,n,1);

        objfun = 2 + sum(((x'.*indexvarN)*indexvarMt-1).^2);

        lb = -inf;
        ub = inf;


    case 91
        %arwhead    求解失败   12957837  
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);

        objfun = sum(-4*term1+3) + sum(((term1.^2) + x(n)^2 ).^2);

        lb = -inf;
        ub = inf;

    case 92
        %bard       12868094 
        n = 3;
        m = 15;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        u = 1:m;
        v = 16-u;
        w = min(u,v);
        y = [0.14 0.18 0.22 0.25 0.29 0.32 ...
             0.35 0.39 0.37	0.58 0.73 0.96 ...
             1.34 2.10 4.39];
        ut = u';
        vt = v';
        wt = w';
        yt = y';
        x1 = x(1) * ones(m,1);
        x2 = x(2) * ones(m,1);
        x3 = x(3) * ones(m,1);

        objfun = sum((yt - (x1 + ut./(vt.*x2+wt.*x3))).^2);

        lb = -inf;
        ub = inf;

    case 93
        %bdexp      12957846    内存不足
        n = 1000;
        ngs = n-2;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:ngs);
        term2 = x(2:ngs+1);
        term3 = x(3:ngs+2);

        objfun = sum( (term1 + term2) ...
                    .* exp((term1+term2).*(-term3)) );


        lb = -inf;
        ub = inf;


    case 94
        %bdqrtic        12868098
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);
        
        term1 = x(1:n-4);
        term2 = x(2:n-3);
        term3 = x(3:n-2);
        term4 = x(4:n-1);

        objfun = sum( (-4*term1+3).^2 ) ...
                + sum( (((term1.^2) ...
                        + 2*(term2.^2) + 3*(term3.^2) ...
                        + 4*(term4.^2) + 5*x(n)^2).^2) );
       

        lb = -inf;
        ub = inf;

    case 95
        %biggs6  对了     12868107
        n = 6;
        m = 13;
        x0 = [1 2 1 1 4 3]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1)*ones(m,1);
        x2 = x(2)*ones(m,1);
        x3 = x(3)*ones(m,1);
        x4 = x(4)*ones(m,1);
        x5 = x(5)*ones(m,1);
        x6 = x(6)*ones(m,1);

        objfun = sum( ( (-exp(-0.1*ones(m,1))) ...
                    +(5*exp(-ones(m,1))) ...
                    -(3*exp(-0.4*ones(m,1))) ...
                    + x3.*exp(-0.1*ones(m,1).*x1) ...
                    - x4.*exp(-0.1*ones(m,1).*x2) ...
                    + x6.*exp(-0.1*ones(m,1).*x5)).^2 );


        lb = -inf;
        ub = inf;

    case 96
        %box3   正确  12868149
        n = 3;
        m = 10;
        x0 = [0 10 1]';

        % Define the objective function
        x = optimvar('x',n,1);

        t = 0.1*(1:m)';
        x1 = x(1)*ones(m,1);
        x2 = x(2)*ones(m,1);
        x3 = x(3)*ones(m,1);

        objfun = sum( (exp(-t.*x1) - exp(-t.*x2) ...
                        -x3.*exp(-t) ...
                        +x3.*exp(-10*t)).^2 );

        lb = -inf;
        ub = inf;

    case 97
        %brkmcc     正确 12868177
        n = 2;
        x0 = [2 2]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = (x1-2)^2 + (x2-1)^2 ... 
                + (1/(1-0.25*x1^2-x2^2))/25 ...
                + 5*(x1-2*x2+1)^2;

        lb = -inf;
        ub = inf;

    case 98
        %brownal    正确  12868190
        n = 10;
        x0 = 1/2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);
        term2 = sum(x);
        term3 = (prod(x)-1)^2;

        objfun = sum( (term1 + term2*ones(9,1) - 11).^2 ) + term3;

        lb = -inf;
        ub = inf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 99
        %brownbs    CNMTr求不出初始点 12872617
        n = 2;
        x0 = ones(n,1);
%         x0 = zeros(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

%         term1 = x(1:n-1);
%         term2 = x(2:n);
% 
%         objfun = sum((term1 - 1000000.0).^2) ...
%             + sum((term2 - 0.000002).^2) ...
%             + sum((term1.*term2 - 2.0).^2);

        x1 = x(1);
        x2 = x(2);

        objfun = (x1 - 1000000.0)^2 ...
                +(x2 - 0.000002)^2 ...
                +(x1*x2 - 2.0)^2;

        lb = -inf;
        ub = inf;


    case 100
        %brownden  12872716  
        n = 4;
        m = 20;
        x0 = [25 5 -5 -1]';

        % Define the objective function
        x = optimvar('x',n,1);
        t = 1:20;
        tt = (t/5)';
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);

        objfun = sum( ((x1+tt*x2-exp(tt)).^2 + (x3+x4*sin(tt)-cos(tt)).^2).^2 );

        lb = -inf;
        ub = inf;

    case 101
        %broydn7d  12872726
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(2:n-1);
        term2 = x(1:n-2);
        term3 = x(3:n);
        term4 = x(1:(n/2));
        term5 = x((n/2+1):n);

        termVer1 = (-2*x(2)+1+(3-2*x(1))*x(1));
        termVer2 = (1-term2 - 2*term3 + (3-2*term1).*term1);
        termVer3 = (-x(n-1)+1+ (3-2*x(n)) *x(n));
        termVer4 = (term4 + term5);


        objfun = (sqrt(termVer1^2))^(7/3) ...
            + sum( (sqrt(termVer2.^2)).^(7/3) )...
            + (sqrt(termVer3^2))^(7/3)...
            + sum( (sqrt(termVer4.^2)).^(7/3) );

        lb = -inf;
        ub = inf;

    case 102
        %chainwoo   12872931
        ns = 499;
        n = 2*ns + 2;
        x0 = -2*ones(n,1);
        x0(1) = -3;
        x0(2) = -1;
        x0(3) = -3;
        x0(4) = -1;

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:2:n-2); %x[2*i-1]
        term2 = x(2:2:n-2); %x[2*i]
        term3 = x(3:2:n);   %x[2*i+1]
        term4 = x(4:2:n);   %x[2*i+2]

        objfun = 1 + ...
        sum ( 100*((term2-term1.^2).^2) + ...
            ((1.0-term1).^2) + ...
            90*((term4-term3.^2).^2) + ...
            ((1.0-term3).^2) + ...
            10*((term2+term4-2.0).^2) + ...
            (((term2-term4).^2)./10) );

        lb = -inf;
        ub = inf;

       
    case 103
        %chnrosnb   12878663
        n = 50;
        x0 = -ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);
        alph = [1.25 1.40 2.40 1.40 1.75 1.20 ...
                2.25 1.20 1.00 1.10 1.50 1.60 ...
                1.25 1.25 1.20 1.20 1.40 0.50 ...
                0.50 1.25 1.80 0.75 1.25 1.40 ...
                1.60 2.00 1.00 1.60 1.25 2.75 ...
                1.25 1.25 1.25 3.00 1.50 2.00 ...
                1.25 1.40 1.80 1.50 2.20 1.40 ...
                1.50 1.25 2.00 1.50 1.25 1.40 ...
                0.60 1.50]';

        term1 = x(1:n-1);
        term2 = x(2:n);
        term3 = alph(2:n);

        objfun = sum (...
    	((term1-(term2.^2)).^2)*16.*(term3.^2) + ...
    	((term2-1.0).^2) );

        lb = -inf;
        ub = inf;

    case 104
        %cosine  12878676       
        n = 1000;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);
        term2 = x(2:n);

        objfun = sum( cos(-0.5*term2+(term1.^2)) );

        lb = -inf;
        ub = inf;

    case 105
        %cragglvy  12878765(1000) 12878694(5000)
        m = 499;
        n = 2*m+2;
        x0 = 2*ones(n,1);
        x0(1) = 1;

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:2:n-2); %2i-1
        term2 = x(2:2:n-2); %2i
        term3 = x(3:2:n);   %2i+1
        term4 = x(4:2:n);   %2i+2

        objfun = sum(...
    	((exp(term1)-term2).^4) + ...
    	100*((term2-term3).^6) + ...
    	((tan(term3-term4)+term3-term4).^4) + ...
    	((term1).^8) + ...
    	((term4-1.0).^2) );

        lb = -inf;
        ub = inf;

    case 106
        %cube   12880725
        n = 2;
        x0 = [-1.2 1]';

        % Define the objective function
        x = optimvar('x',n,1);
        x1 = x(1);
        x2 = x(2);
        term1 = x(2:n);
        term2 = x(1:n-1);

        objfun = (x1-1)^2 + ...
            sum( 100* ((term1- (term2).^3).^2) );
       
%         objfun = (x1-1)^2 + ...
%                  100*(x2-x1^3)^2;


        lb = -inf;
        ub = inf;

%     case 107
%         %curly10  12880844    %自动微分不支持变量作为指数
%         n = 1000;
%         k = 10;
%         xi = 1:n;
%         x0 = (0.0001*xi/(n+1))';
% 
%         % Define the objective function
%         x = optimvar('x',n,1);
% 
%         X = repmat(x,1,n);
% %         Xtr = tril(X);  %下三角 不能用该函数
%         e = tril(ones(n));
%         Xtr = X.*e;
%         XtrNk = [Xtr(:,11:n),zeros(n,10)];
% 
%         objfun = sum(sum(Xtr-XtrNk));
% 
%         lb = -inf;
%         ub = inf;

    case 108
        %dqdrtic    12882615
        n = 1000;
        x0 = 3*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-2);   %xi
        term2 = x(2:n-1);   %xi+1
        term3 = x(3:n);     %xi+2   

        objfun = sum((100*(term2.^2)+100*(term3.^2)+(term1.^2)));

        lb = -inf;
        ub = inf;

    case 109
        %dqrtic  12882616
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        var = (1:n)';

        objfun = sum( (x-var).^4 );

        lb = -inf;
        ub = inf;

    case 110
        %eg2        12882622
        n = 1000;
        x0 = zeros(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);   %xi

        objfun = sum( sin(x(1) + (term1.^2) - 1.0) )...
            + 0.5*sin(x(n)^2);

        lb = -inf;
        ub = inf;

    case 111
        %engval1     12882629  
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);   %xi
        term2 = x(2:n);     %xi+1

        objfun = sum( ((term1.^2)+(term2.^2)).^2 ) ...
                + sum( (-4*term1+3.0) );

        lb = -inf;
        ub = inf;

    case 112
        %engval2        12882702
        n = 3;
        x0 = [1 2 0]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        x3 = x(3);

        objfun = (x1^2+x2^2+x3^2-1)^2   ...
    	+ (x1^2+x2^2+(x3-2)^2-1)^2  ...
    	+ (x1+x2+x3-1)^2    ...
    	+ (x1+x2-x3+1)^2    ....
    	+ (3*x2^2+x1^3+(5*x3-x1+1)^2-36)^2;

        lb = -inf;
        ub = inf;

    case 113
        %errinros   12882712
        n = 50;
        x0 = -ones(n,1);
        alpha = [1.25 1.40 2.40 1.40 1.75 1.20 ...
                 2.25 1.20 1.00 1.10 1.50 1.60 ...
                 1.25 1.25 1.20 1.20 1.40 0.50 ...
                 0.50 1.25 1.80 0.75 1.25 1.40 ...
                 1.60 2.00 1.00 1.60 1.25 2.75 ...
                 1.25 1.25 1.25 3.00 1.50 2.00 ...
                 1.25 1.40 1.80 1.50 2.20 1.40 ...
                 1.50 1.25 2.00 1.50 1.25 1.40 ...
                 0.60 1.50]';


        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);   %xi-1
        term2 = x(2:n);     %xi
        term3 = alpha(2:n);

        objfun = sum( ( term1-((16*(term3.^2)).*(term2.^2)) ).^2 ) ...
        	+sum( (term2-1.0).^2 );

        lb = -inf;
        ub = inf;

    case 114
        %extrosnb  12883028
        n = 10;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(2:n);   %xi
        term2 = x(1:n-1);     %xi+1

        objfun = (x(1)-1)^2 + ...
                sum( 100* ((term1-(term2.^2)).^2) );

        lb = -inf;
        ub = inf;

    case 115
        %fletcbv3  12883031
        n = 1000;
        var = (1:n)';
        kappa = 1;
        objscale = 1e+8;
        h = 1/(n+1);
        p = 1/objscale;
        x0 = var/(h);

        % Define the objective function
        x = optimvar('x',n,1);
        

        term1 = x(2:n);   %xi
        term2 = x(1:n-1);    

        objfun = 0.5*p*(x(1)^2 + ...
                sum( 0.5*p*( (term2-term1).^2 ) ) + 0.5*p*x(n)^2  +...
                sum( p*(-1 - 2/h^2)*x ) + ...
                sum( -kappa*p*cos(x)/h^2 )...
                );

        lb = -inf;
        ub = inf;

    case 116
        %genhumps   12883049 
        n = 5;
        zeta = 2;
        x0 = 506.2*ones(n,1);
        x0(1) = -506;

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);   %xi
        term2 = x(2:n);     %xi+1

        objfun = sum( (((sin(zeta*term1)).^2).*((sin(zeta*term2)).^2)) ...
                +0.05*((term1.^2)+(term2.^2)) );

        lb = -inf;
        ub = inf;

%     case 117      %自动微分不支持变量作为指数
%         %gulf       12883074
%         n = 3;
%         m = 99;
%         x0 = [5 2.5 0.15]';
% 
%         % Define the objective function
%         x = optimvar('x',n,1);
% 
%         var = (1:m)';
%         t = var/100;
%         y = 25+ ( (-50*log(t)).^(2/3) );
%         sqr1 = (y-x(2)*ones(m,1)).^2;
%         abs1 = sqrt(sqr1);
%         x3 = x(3);
%         exp3 = power((exp( abs1 )), x3);
% 
% 
%         objfun = sum( (((exp3) / (-x(1)) -t).^2) );
% 
%         lb = -inf;
%         ub = inf;

    case 118
        %hilbertb   12886434
        n = 50;
        d = 5;
        x0 = -3*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        
        e = triu(ones(n));
        mar1 = (1:n)'.*e;           %[1; 1 2; 1 2 3]    列向量.*上三角
        mar2 = repmat((1:n),n,1);   %[1 1 1; 2 2 2; 3 3 3]  %行向量
        mar3 = mar2 .* e;   %1 2 3上三角

        Xtriu0 = repmat(x',n,1);
        Xtriu1 = Xtriu0.*e;              %[x1; x2 x2; x3 x3 x3]    

        Xtriu2 = repmat(x,1,n);     %[x1 x2 x3; x1 x2 x3; x1 x2 x3] 列向量
        xij = (Xtriu1.*Xtriu2)./(mar1+mar3-1);
        sumXij = sum(xij);

        term = (1:n)';

        mul = d+(1./(4*term-2));
        sum2 = (x.^2).*mul;

        objfun = sum( sumXij' + sum2 );

        lb = -inf;
        ub = inf;


    case 119
        %himmelbh  12886485
        n = 2;
        x0 = [0 2]';

        % Define the objective function
        x = optimvar('x',n,1);


        objfun = -3*x(1)-2*x(2)+2+x(1)^3+x(2)^2;

        lb = -inf;
        ub = inf;

    case 120
        %kowosb  12886632
        n = 4;
        m = 11;
        x0 = [0.25 0.39 0.415 0.39]';

        % Define the objective function
        x = optimvar('x',n,1);

        y = [0.1957	0.1947	0.1735	0.1600	0.0844 0.0627 ...
            0.0456	0.0342	0.0323	0.0235	0.0246]';
        u = [4.0000	2.0000	1.0000	0.5000	0.2500 0.1670 ...
            0.1250	0.1000	0.0833	0.0714	0.0625]';

        x1 = x(1)*ones(m,1);
        x2 = x(2)*ones(m,1);
        x3 = x(3)*ones(m,1);
        x4 = x(4)*ones(m,1);


        objfun = sum( ( y-((x1.*((u.^2)+(u.*x2))) ...
                    ./((u.^2)+(u.*x3)+x4)) ).^2 );


        lb = -inf;
        ub = inf;

  %%%%%%%%%%%%%%%%%%%%%%%
    case 121    %再看看
        %meyer3    12921903
        n = 3;
        m = 16;     
        x0 = [0.02 4000 250]';

        % Define the objective function
        x = optimvar('x',n,1);

        t = 45 + 5*(1:m)';
        y = [34780 28610 23650 19630 ...
             16370 13720 11540 9744  ...
             8261  7030  6005  5147  ...
             4427  3820  3307  2872]';
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);

        objfun = sum( (x1*exp(x2./(t+x3)) ...
            -y).^2 );


        lb = -inf;
        ub = inf;

    case 122
        %noncvxu2    12921461
        n = 1000;
        var = (1:n)';
        x0 = var;

        % Define the objective function
        x = optimvar('x',n,1);

        var1 = mod(3*var-2,1000)+1;
        var2 = mod(7*var-3,1000)+1;

        term1 = x(var1);
        term2 = x(var2);
    
        objfun = sum( (x + term1 + term2).^2 ...
            + 4*cos(x + term1 + term2) );


        lb = -inf;
        ub = inf;

    case 123
        %osbornea  12921717
        n = 5;
        m = 33;
        var = (1:m)';
        x0 = [0.5 1.5 -1 0.01 0.02]';

        % Define the objective function
        x = optimvar('x',n,1);

        t = 10*(var-1);
        y = [0.844 0.908 0.932 0.936 ...
             0.925 0.908 0.881 0.850 ...
             0.818 0.784 0.751 0.718 ...
             0.685 0.658 0.628 0.603 ...
             0.580 0.558 0.538 0.522 ...
             0.506 0.490 0.478 0.467 ...
             0.457 0.448 0.438 0.431 ...
             0.424 0.420 0.414 0.411 0.406]';
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        x5 = x(5);

        objfun = sum( ( y-x1-x2*exp(-t*x4) ...
            -x3*exp(-t*x5) ).^2 );


        lb = -inf;
        ub = inf;

    case 124
        %penalty2   12921789
        n = 100;
        m = 2*n;
        x0 = 1/2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        a = 10^-5;
        var = (1:n)';
        y = exp(var/10) + exp((var-1)/10);
        termY = y(2:n);
        termX1 = x(2:n);
        termX2 = x(1:n-1);

        objfun = (x(1)-0.2)^2 + ...
            sum( a*((exp(termX1/10) + exp(termX2/10)-termY).^2) ) + ...
            sum( a*((exp(termX1/10)-exp(-1/10)).^2) ) + ...
            ( sum(  (n-var+1).*(x.^2) )-1 )^2;


        lb = -inf;
        ub = inf;

    case 125
        %quartc  12921867
        n = 1000;
        x0 = 2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        var = (1:n)';

        objfun = sum( (x-var).^4 );

        lb = -inf;
        ub = inf;

    case 126
        %rosenbr  12922047
        n = 2;
        x0 = [-1.2 1.0]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = (x2-x1^2)^2/0.01+(x1-1)^2;

        lb = -inf;
        ub = inf;

    case 127
        %scosine  12922352
        n = 1000;
        var = (1:n)';
        scal = 12.0;
        scale = exp((var-1)*scal/(n-1));
        x0 = 1.0./scale;

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:n-1);
        term2 = x(2:n);
        termS1 = scale(1:n-1);
        termS2 = scale(2:n);

        objfun = sum( cos(-0.5*(termS2.*term2) ...
            + ((termS1.^2).*(term1.^2))) );

        lb = -inf;
        ub = inf;

    case 128
        %scurly10   12928769 原为1w，1000的内存也不够
        n = 100;    %12933177
        k = 10;
        sc = 12;
        var1 = (0:n-1)';
        var2 = (1:n)';
        scale = exp(var1*sc/(n-1));
        x0 = 0.0001*(var2.*scale)/(n+1);

        % Define the objective function
        x = optimvar('x',n,1);

        X = repmat(x,1,n);  
%         Xtr = tril(X);  %下三角 不能用该函数
        e = tril(ones(n));
        Xtr = X.*e;
        Mar0 = Xtr(:,12:n);
        Mar0 = [Mar0,zeros(n,1)];
        Mar1 = Xtr(:,1:n-k)-Mar0;    %下三角的前n-10个
        Mar2 = Xtr(:,n-k+1:n);       %下三角的后10个

        term1 = scale(1:n-k);
        term2 = scale(n-k+1:n);

        sumMar1 = sum(Mar1)';
        y = sumMar1.*term1;
        z = (sum(Mar2))'.*term2;

        objfun = sum( y.*((y.*((y.^2)-20))-0.1) ) +...
                 sum( z.*((z.*((z.^2)-20))-0.1) );

        lb = -inf;
        ub = inf;    

    case 129
        %sineval    12933707
        n = 2;
        x0 = [4.712389 -1.0]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);
        c = 1e-4;

        objfun = (x2-sin(x1))^2/c + x1^2/4;

        lb = -inf;
        ub = inf;   


    case 130
        %sinquad 12933733
        n = 1000;
        x0 = 0.1*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        term = x(2:n-1);
        x1 = x(1);
        xn = x(n);


        objfun = (x1-1)^4 ...
        	+ sum( (sin(term-xn)-x1^2+(term.^2)).^2 ) ...
        	+ (xn^2-x1^2)^2	;

        lb = -inf;
        ub = inf;   

    case 131
        %sisser     12933761
        n = 2;
        x0 = [1.0 0.1]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = 3*x1^4 - ...
                 2*(x1*x2)^2 + 3*x2^4;

        lb = -inf;
        ub = inf;    

    case 132
        %srosenbr     12933776
        n = 1000;
        x0 = ones(n,1);
        x0(1:2:n) = -1.2*x0(1:2:n);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(2:2:n);   %x2i
        term2 = x(1:2:n);   %x2i-1

        objfun = sum( 100*((term1-(term2.^2)).^2) ...
                    + ((term2-1).^2) );

        lb = -inf;
        ub = inf; 

    case 133
        %tridia     12933818    
        n = 1000;       %problem.objective算的不一样；
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        alpha = 2.0;
        beta = 1.0;
        gamma = 1.0;
        delta = 1.0;

        term1 = x(1:n-1);   %xi-1
        term2 = x(2:n);   %xi
        var = (2:n)';

        objfun = gamma*(x(1)*delta-1.0)^2 ...
                + sum( var.*((-beta*term1+alpha*term2).^2) );

        lb = -inf;
        ub = inf; 

    case 134
        %vardim     12934785 
        n = 100;
        var = (1:n)';
        x0 = 1-(var/n);

        % Define the objective function
        x = optimvar('x',n,1);

        objfun = sum( (x-1).^2 ) + ...
                 (sum(var.*x) - n*(n+1)/2)^2 + ...
                 (sum(var.*x) - n*(n+1)/2)^4;

        lb = -inf;
        ub = inf;  

    case 135
        %watson     12934816
        n = 31;
        x0 = zeros(31,1);

        % Define the objective function
        x = optimvar('x',n,1);

        var1 = (1:29)';
        t = var1/29;

        var2 = (1:31)';
        term1 = var2(1:30); %term1=1:30 j-1

        x1 = x(1);
        x2 = x(2);

        mar1 = repmat(t',n,1);      %[t1 t2...t29; t1 t2...t29...]  31行29列(1~31)
        mar2 = repmat(t',n-1,1);    %[t1 t2...t29; t1 t2...t29...]  30行29列(2~31)
        term5 = x(2:31);       %行向量的x，做向量*矩阵
        
        v1 = (0:29)';    %j-2
        mv1 = repmat(v1,1,29);
        v2 = (0:30)';    %j-1
        mv2 = repmat(v2,1,29);

        part1 = ((term1.*term5)')*(mar2.^mv1);
        part2 = x'*(mar1.^mv2);

%         objfun = sum {i in 1..29}
%         (sum {j in 2..31} (j-1)*x[j]*t[i]^(j-2) ...
%             - (sum {j in 1..31} x[j]*t[i]^(j-1))^2 - 1)^2 ... 
        objfun = sum( (part1 - part2.^2 -1).^2) ... 
        + x1^2 + (x2-x1^2-1)^2;

        lb = -inf;
        ub = inf;  

     case 136
        %woods  12939123
        n = 1000;
        x0 = -ones(n,1);
        x0(1:2:n) = -3*x0(1:2:n);

        % Define the objective function
        x = optimvar('x',n,1);

        term1 = x(1:4:n);   %x4i-3
        term2 = x(2:4:n);   %x4i-2
        term3 = x(3:4:n);   %x4i-1
        term4 = x(4:4:n);   %x4i

        objfun = sum( 100*((term2-(term1.^2)).^2) + ...
        ((1-term1).^2) + ...
        90*((term4-(term3.^2)).^2) + ...
        ((1-term3).^2) + ...
        10*((term2+term4-2).^2) + ...
        0.1*((term2-term4).^2) );


        lb = -inf;
        ub = inf;    

     case 137
        %zangwil2   12939147
        n = 2;
        x0 = [3 8]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = (-56*x1-256*x2+991+ ...
            16*x1^2+16*x2^2- ...
            8*x1*x2)/15;

        lb = -inf;
        ub = inf;  

    case 138
        %fletchbv   12939315
        n = 1000;
        kappa = 1.0;
        objscale = 1.0D+0;
        h = 1/(n+1);
        p = 1/objscale;
        var = (1:n)';
        x0 = var*h;

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        xn = x(n);

        term1 = x(1:n-1);   %x1 - xn-1
        term2 = x(2:n);

        objfun = 0.5*p*(x1)^2 + ...
    	sum( 0.5*p*((term1-term2).^2) )  +  ...
    	0.5*p*(xn)^2 +  ...
    	sum  ( p*(-1-2/h^2)*x ) +   ...
    	sum  ( -kappa*p*cos(x)/h^2 );

        lb = -inf;
        ub = inf; 

    case 139
        %heart6ls   12945434
        n = 6;
        x0 = [0 0 1 1 1 1]';

        % Define the objective function
        x = optimvar('x',n,1);

        a = x(1);
        c = x(2);
        t = x(3);
        u = x(4);
        v = x(5);
        w = x(6);

        objfun = (t * a + u * (-0.816-a) - v * c - w * (-0.017-c) + 1.826)^2 + (v * a + w * ...
        	(-0.816-a) + t * c + u * (-0.017-c) + 0.754)^2 + (a * (t^2-v^2) - 2.0*c * t * v ...
        	+ (-0.816-a) * (u^2-w^2) - 2.0*(-0.017-c) * u * w + 4.839)^2 + (c * (t^2-v^2) + ...
        	2.0*a * t * v + (-0.017-c) * (u^2-w^2) + 2.0*(-0.816-a) * u * w + 3.259)^2 + (a ...
        	* t * (t^2-(3.0)*v^2) + c * v * (v^2-(3.0)*t^2) + (-0.816-a) * u *              ...
        	(u^2-(3.0)*w^2) + (-0.017-c) * w * (w^2-(3.0)*u^2) + 14.023)^2 + (c * t *       ...
        	(t^2-(3.0)*v^2) - a * v * (v^2-(3.0)*t^2) + (-0.017-c) * u * (u^2-(3.0)*w^2) -  ...
        	(-0.816-a) * w * (w^2-(3.0)*u^2) - 15.467)^2;

        lb = -inf;
        ub = inf;

    case 140
        %heart8ls   12947616
        n = 8;
        x0 = ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        a = x(1);
        b = x(2);
        c = x(3);
        d = x(4);
        t = x(5);
        u = x(6);
        v = x(7);
        w = x(8);

        objfun = (a + b + 0.69)*(a + b + 0.69) + (c + d + 0.044)*(c + d + 0.044) + (t * a + u * ...
        	b - v * c - w * d + 1.57)*(t * a + u * b - v * c - w * d + 1.57) + (v * a + w * ...
        	b + t * c + u * d + 1.31)*(v * a + w * b + t * c + u * d + 1.31) + (a * ...
        	(t^2-v^2) - 2.0*c * t * v + b * (u^2-w^2) - 2.0*d * u * w + 2.65)*(a * ...
        	(t^2-v^2) - 2.0*c * t * v + b * (u^2-w^2) - 2.0*d * u * w + 2.65) + (c * ...
        	(t^2-v^2) + 2.0*a * t * v + d * (u^2-w^2) + 2.0*b * u * w - 2.0)*(c * (t^2-v^2) ...
        	+ 2.0*a * t * v + d * (u^2-w^2) + 2.0*b * u * w - 2.0) + (a * t * ...
        	(t^2-(3.0)*v^2) + c * v * (v^2-(3.0)*t^2) + b * u * (u^2-(3.0)*w^2) + d * w * ...
        	(w^2-(3.0)*u^2) + 12.6)*(a * t * (t^2-(3.0)*v^2) + c * v * (v^2-(3.0)*t^2) + b ...
        	* u * (u^2-(3.0)*w^2) + d * w * (w^2-(3.0)*u^2) + 12.6) + (c * t * ...
        	(t^2-(3.0)*v^2) - a * v * (v^2-(3.0)*t^2) + d * u * (u^2-(3.0)*w^2) - b * w * ...
        	(w^2-(3.0)*u^2) - 9.48)*(c * t * (t^2-(3.0)*v^2) - a * v * (v^2-(3.0)*t^2) + d ...
        	* u * (u^2-(3.0)*w^2) - b * w * (w^2-(3.0)*u^2) - 9.48);
        lb = -inf;
        ub = inf;

    case 141
        %hilberta   12947648
        n = 10;
        x0 = ones(n,1);
        x0(1) = -4;
        x0(2) = -2;

        % Define the objective function
        x = optimvar('x',n,1);

        A = zeros(n,n);

        for i=1:n
            for j=1:n
                A(i,j) = 1/(i+j-1);
            end
        end

        objfun = x'*A*x;

        lb = -inf;
        ub = inf;  

    case 142
        %himmelbf   2947679
        n = 4;
        x0 = [2.7 90 1500 10]';

        % Define the objective function
        x = optimvar('x',n,1);

        a = [0.0 0.000428 0.001000 0.001610 ...
            0.002090 0.003480 0.005250]';
        b = [7.391 11.18 16.44 16.20 22.20 ...
             24.02 31.32]';

        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);

        objfun = 10000*sum( ...
            (-1 + ...
            ((x1^2+a*x2^2+(a.^2)*x3^2) ./ ...
            (b.*(1+a*x4^2))) ...
            ).^2 ...
            );

        lb = -inf;
        ub = inf;  

    case 143
        %humps  12947711
        n = 2;
        x0 = [-506.0 -506.2]';

        % Define the objective function
        x = optimvar('x',n,1);

        zeta = 20.0;
        x1 = x(1);
        y = x(2);

        objfun = 0.05*(x1^2+y^2) + ...
            (sin(zeta*x1)*sin(zeta*y))^2;

        lb = -inf;
        ub = inf; 

    case 144
        %methanb8   12947753 
        n = 31;
        x0 = [107.47 0.09203 0.908 102.4 0.1819 ...
              0.8181 97.44 0.284 0.716 96.3 ...
              0.3051 0.6949 93.99 0.3566 0.6434 ...
              89.72 0.468 0.532 83.71 0.6579 ...
              0.3421 78.31 0.8763 0.1237 886.37 ...
              910.01 922.52 926.46 935.56 952.83 ...
              975.73]';

        % Define the objective function
        x = optimvar('x',n,1);

        t0  = x(1);
        x0_1 = x(2);
    	x0_2 = x(3);
    	t1 = x(4);
    	x1_1 = x(5);
    	x1_2 = x(6);
    	t2 = x(7);
    	x2_1 = x(8);
    	x2_2 = x(9);
    	t3 = x(10);
    	x3_1 = x(11);
    	x3_2 = x(12);
    	t4 = x(13);
    	x4_1 = x(14);
    	x4_2 = x(15);
    	t5 = x(16);
    	x5_1 = x(17);
    	x5_2 = x(18);
    	t6 = x(19);
    	x6_1 = x(20);
    	x6_2 = x(21);
    	t7 = x(22);
    	x7_1 = x(23);
    	x7_2 = x(24);
    	v0 = x(25);
    	v1 = x(26);
    	v2 = x(27);
    	v3 = x(28);
    	v4 = x(29);
    	v5 = x(30);
    	v6 = x(31);


        objfun = (((-1.0 * x1_1 * (v0 + 693.37 )  + ...
	(v0*x0_1*1.0*8.264462809917355e-4*exp(18.5751+(-3632.649/(t0+239.2)))) + ...
	693.37*x0_1)/100.0))^2 + ...
	((x6_1*1.0*8.695652173913045e-4*exp(18.5751+(-3632.649/(t6+239.2)))) - x7_1)^2 ...
	+ (((-1.0 * x2_1 * (v1 + 693.37 )  + ...
	(v0*x0_1*-1.0*8.264462809917355e-4*exp(18.5751+(-3632.649/(t0+239.2)))) + 1.0 * ...
	x1_1 * (v0 + 693.37 )  + ...
	(v1*x1_1*1.0*8.333333333333334e-4*exp(18.5751+(-3632.649/(t1+239.2)))))/100.0))^2 + (((-1.0 * x3_1 * (v2 + -442.13 )  + ...
	(v1*x1_1*-1.0*8.333333333333334e-4*exp(18.5751+(-3632.649/(t1+239.2)))) + 1.0 * ...
	x2_1 * (v1 + 693.37 )  + ...
	(v2*x2_1*1.0*8.403361344537816e-4*exp(18.5751+(-3632.649/(t2+239.2)))) - ...
	451.25)/100.0))^2 + (((-1.0 * x4_1 * (v3 + -442.13 )  + ...
	(v2*x2_1*-1.0*8.403361344537816e-4*exp(18.5751+(-3632.649/(t2+239.2)))) + 1.0 * ...
	x3_1 * (v2 + -442.13 )  + ...
	(v3*x3_1*1.0*8.474576271186439e-4*exp(18.5751+(-3632.649/(t3+239.2)))))/100.0))^2 + (((-1.0 * x5_1 * (v4 + -442.13 )  + ...
	(v3*x3_1*-1.0*8.474576271186439e-4*exp(18.5751+(-3632.649/(t3+239.2)))) + 1.0 * ...
	x4_1 * (v3 + -442.13 )  + ...
	(v4*x4_1*1.0*8.547008547008547e-4*exp(18.5751+(-3632.649/(t4+239.2)))))/100.0))^2 + (((-1.0 * x6_1 * (v5 + -442.13 )  + ...
	(v4*x4_1*-1.0*8.547008547008547e-4*exp(18.5751+(-3632.649/(t4+239.2)))) + 1.0 * ...
	x5_1 * (v4 + -442.13 )  + ...
	(v5*x5_1*1.0*8.620689655172415e-4*exp(18.5751+(-3632.649/(t5+239.2)))))/100.0))^2 + (((-1.0 * x7_1 * (v6 + -442.13 )  + ...
	(v5*x5_1*-1.0*8.620689655172415e-4*exp(18.5751+(-3632.649/(t5+239.2)))) + 1.0 * ...
	x6_1 * (v5 + -442.13 )  + ...
	(v6*x6_1*1.0*8.695652173913045e-4*exp(18.5751+(-3632.649/(t6+239.2)))))/100.0))^2 + (((-1.0 * x1_2 * (v0 + 693.37 )  + ...
	(v0*x0_2*1.0*8.264462809917355e-4*exp(18.3443+(-3841.2203/(t0+228.0)))) + ...
	693.37*x0_2)/100.0))^2 + ...
	((x6_2*1.0*8.695652173913045e-4*exp(18.3443+(-3841.2203/(t6+228.0)))) - x7_2)^2 ...
	+ (((-1.0 * x2_2 * (v1 + 693.37 )  + ...
	(v0*x0_2*-1.0*8.264462809917355e-4*exp(18.3443+(-3841.2203/(t0+228.0)))) + 1.0 ...
	* x1_2 * (v0 + 693.37 )  + ...
	(v1*x1_2*1.0*8.333333333333334e-4*exp(18.3443+(-3841.2203/(t1+228.0)))))/100.0))^2 + (((-1.0 * x3_2 * (v2 + -442.13 )  + ...
	(v1*x1_2*-1.0*8.333333333333334e-4*exp(18.3443+(-3841.2203/(t1+228.0)))) + 1.0 ...
	* x2_2 * (v1 + 693.37 )  + ...
	(v2*x2_2*1.0*8.403361344537816e-4*exp(18.3443+(-3841.2203/(t2+228.0)))) - ...
	684.25)/100.0))^2 + (((-1.0 * x4_2 * (v3 + -442.13 )  + ...
	(v2*x2_2*-1.0*8.403361344537816e-4*exp(18.3443+(-3841.2203/(t2+228.0)))) + 1.0 ...
	* x3_2 * (v2 + -442.13 )  + ...
	(v3*x3_2*1.0*8.474576271186439e-4*exp(18.3443+(-3841.2203/(t3+228.0)))))/100.0))^2 + (((-1.0 * x5_2 * (v4 + -442.13 )  + ...
	(v3*x3_2*-1.0*8.474576271186439e-4*exp(18.3443+(-3841.2203/(t3+228.0)))) + 1.0 ...
	* x4_2 * (v3 + -442.13 )  + ...
	(v4*x4_2*1.0*8.547008547008547e-4*exp(18.3443+(-3841.2203/(t4+228.0)))))/100.0))^2 + (((-1.0 * x6_2 * (v5 + -442.13 )  + ...
	(v4*x4_2*-1.0*8.547008547008547e-4*exp(18.3443+(-3841.2203/(t4+228.0)))) + 1.0 ...
	* x5_2 * (v4 + -442.13 )  + ...
	(v5*x5_2*1.0*8.620689655172415e-4*exp(18.3443+(-3841.2203/(t5+228.0)))))/100.0))^2 + (((-1.0 * x7_2 * (v6 + -442.13 )  + ...
	(v5*x5_2*-1.0*8.620689655172415e-4*exp(18.3443+(-3841.2203/(t5+228.0)))) + 1.0 ...
	* x6_2 * (v5 + -442.13 )  + ...
	(v6*x6_2*1.0*8.695652173913045e-4*exp(18.3443+(-3841.2203/(t6+228.0)))))/100.0))^2 + ((x0_1*1.0*8.264462809917355e-4*exp(18.5751+(-3632.649/(t0+239.2)))) + ...
	(x0_2*1.0*8.264462809917355e-4*exp(18.3443+(-3841.2203/(t0+228.0)))) - 1.0)^2 + ...
	((x1_1*1.0*8.333333333333334e-4*exp(18.5751+(-3632.649/(t1+239.2)))) + ...
	(x1_2*1.0*8.333333333333334e-4*exp(18.3443+(-3841.2203/(t1+228.0)))) - 1.0)^2 + ...
	((x2_1*1.0*8.403361344537816e-4*exp(18.5751+(-3632.649/(t2+239.2)))) + ...
	(x2_2*1.0*8.403361344537816e-4*exp(18.3443+(-3841.2203/(t2+228.0)))) - 1.0)^2 + ...
	((x3_1*1.0*8.474576271186439e-4*exp(18.5751+(-3632.649/(t3+239.2)))) + ...
	(x3_2*1.0*8.474576271186439e-4*exp(18.3443+(-3841.2203/(t3+228.0)))) - 1.0)^2 + ...
	((x4_1*1.0*8.547008547008547e-4*exp(18.5751+(-3632.649/(t4+239.2)))) + ...
	(x4_2*1.0*8.547008547008547e-4*exp(18.3443+(-3841.2203/(t4+228.0)))) - 1.0)^2 + ...
	((x5_1*1.0*8.620689655172415e-4*exp(18.5751+(-3632.649/(t5+239.2)))) + ...
	(x5_2*1.0*8.620689655172415e-4*exp(18.3443+(-3841.2203/(t5+228.0)))) - 1.0)^2 + ...
	((x6_1*1.0*8.695652173913045e-4*exp(18.5751+(-3632.649/(t6+239.2)))) + ...
	(x6_2*1.0*8.695652173913045e-4*exp(18.3443+(-3841.2203/(t6+228.0)))) - 1.0)^2 + ...
	((x7_1*1.0*8.771929824561405e-4*exp(18.5751+(-3632.649/(t7+239.2)))) + ...
	(x7_2*1.0*8.771929824561405e-4*exp(18.3443+(-3841.2203/(t7+228.0)))) - 1.0)^2 + ...
	((((v0*x0_1*1.0*8.264462809917355e-4*exp(18.5751+(-3632.649/(t0+239.2)))) * ...
	(9566.67+-1.59*t0+0.0422*t0*t0) + 693.37 * x0_1 * (0.0+15.97*t0+0.0422*t0*t0) + ...
	-1.0 * x1_1 * (693.37 + v0 ) *(0.0+15.97*t1+0.0422*t1*t1) + ...
	(v0*x0_2*1.0*8.264462809917355e-4*exp(18.3443+(-3841.2203/(t0+228.0)))) * ...
	(10834.67+8.74*t0+0.0*t0*t0) + 693.37 * x0_2 * (0.0+18.1*t0+0.0*t0*t0) + -1.0 * ...
	x1_2 * (693.37 + v0 ) *(0.0+18.1*t1+0.0*t1*t1) - 8386200.0)/100000.0))^2 + ...
	((((v1*x1_1*1.0*8.333333333333334e-4*exp(18.5751+(-3632.649/(t1+239.2)))) * ...
	(9566.67+-1.59*t1+0.0422*t1*t1) + 1.0 * x1_1 * (693.37 + v0 ) ...
	*(0.0+15.97*t1+0.0422*t1*t1) + ...
	(v0*x0_1*-1.0*8.264462809917355e-4*exp(18.5751+(-3632.649/(t0+239.2)))) * ...
	(9566.67+-1.59*t0+0.0422*t0*t0) + -1.0 * x2_1 * (693.37 + v1 ) ...
	*(0.0+15.97*t2+0.0422*t2*t2) + ...
	(v1*x1_2*1.0*8.333333333333334e-4*exp(18.3443+(-3841.2203/(t1+228.0)))) * ...
	(10834.67+8.74*t1+0.0*t1*t1) + 1.0 * x1_2 * (693.37 + v0 ) ...
	*(0.0+18.1*t1+0.0*t1*t1) + ...
	(v0*x0_2*-1.0*8.264462809917355e-4*exp(18.3443+(-3841.2203/(t0+228.0)))) * ...
	(10834.67+8.74*t0+0.0*t0*t0) + -1.0 * x2_2 * (693.37 + v1 ) ...
	*(0.0+18.1*t2+0.0*t2*t2))/100000.0))^2 + ...
	((((v2*x2_1*1.0*8.403361344537816e-4*exp(18.5751+(-3632.649/(t2+239.2)))) * ...
	(9566.67+-1.59*t2+0.0422*t2*t2) + 1.0 * x2_1 * (693.37 + v1 ) ...
	*(0.0+15.97*t2+0.0422*t2*t2) + ...
	(v1*x1_1*-1.0*8.333333333333334e-4*exp(18.5751+(-3632.649/(t1+239.2)))) * ...
	(9566.67+-1.59*t1+0.0422*t1*t1) + -1.0 * x3_1 * (-442.13 + v2 ) ...
	*(0.0+15.97*t3+0.0422*t3*t3) + ...
	(v2*x2_2*1.0*8.403361344537816e-4*exp(18.3443+(-3841.2203/(t2+228.0)))) * ...
	(10834.67+8.74*t2+0.0*t2*t2) + 1.0 * x2_2 * (693.37 + v1 ) ...
	*(0.0+18.1*t2+0.0*t2*t2) + ...
	(v1*x1_2*-1.0*8.333333333333334e-4*exp(18.3443+(-3841.2203/(t1+228.0)))) * ...
	(10834.67+8.74*t1+0.0*t1*t1) + -1.0 * x3_2 * (-442.13 + v2 ) ...
	*(0.0+18.1*t3+0.0*t3*t3) - 1894471.11025)/100000.0))^2 + ...
	((((v3*x3_1*1.0*8.474576271186439e-4*exp(18.5751+(-3632.649/(t3+239.2)))) * ...
	(9566.67+-1.59*t3+0.0422*t3*t3) + 1.0 * x3_1 * (-442.13 + v2 ) ...
	*(0.0+15.97*t3+0.0422*t3*t3) + ...
	(v2*x2_1*-1.0*8.403361344537816e-4*exp(18.5751+(-3632.649/(t2+239.2)))) * ...
	(9566.67+-1.59*t2+0.0422*t2*t2) + -1.0 * x4_1 * (-442.13 + v3 ) ...
	*(0.0+15.97*t4+0.0422*t4*t4) + ...
	(v3*x3_2*1.0*8.474576271186439e-4*exp(18.3443+(-3841.2203/(t3+228.0)))) * ...
	(10834.67+8.74*t3+0.0*t3*t3) + 1.0 * x3_2 * (-442.13 + v2 ) ...
	*(0.0+18.1*t3+0.0*t3*t3) + ...
	(v2*x2_2*-1.0*8.403361344537816e-4*exp(18.3443+(-3841.2203/(t2+228.0)))) * ...
	(10834.67+8.74*t2+0.0*t2*t2) + -1.0 * x4_2 * (-442.13 + v3 ) ...
	*(0.0+18.1*t4+0.0*t4*t4))/100000.0))^2 + ...
	((((v4*x4_1*1.0*8.547008547008547e-4*exp(18.5751+(-3632.649/(t4+239.2)))) * ...
	(9566.67+-1.59*t4+0.0422*t4*t4) + 1.0 * x4_1 * (-442.13 + v3 ) ...
	*(0.0+15.97*t4+0.0422*t4*t4) + ...
	(v3*x3_1*-1.0*8.474576271186439e-4*exp(18.5751+(-3632.649/(t3+239.2)))) * ...
	(9566.67+-1.59*t3+0.0422*t3*t3) + -1.0 * x5_1 * (-442.13 + v4 ) ...
	*(0.0+15.97*t5+0.0422*t5*t5) + ...
	(v4*x4_2*1.0*8.547008547008547e-4*exp(18.3443+(-3841.2203/(t4+228.0)))) * ...
	(10834.67+8.74*t4+0.0*t4*t4) + 1.0 * x4_2 * (-442.13 + v3 ) ...
	*(0.0+18.1*t4+0.0*t4*t4) + ...
	(v3*x3_2*-1.0*8.474576271186439e-4*exp(18.3443+(-3841.2203/(t3+228.0)))) * ...
	(10834.67+8.74*t3+0.0*t3*t3) + -1.0 * x5_2 * (-442.13 + v4 ) ...
	*(0.0+18.1*t5+0.0*t5*t5))/100000.0))^2 + ...
	((((v5*x5_1*1.0*8.620689655172415e-4*exp(18.5751+(-3632.649/(t5+239.2)))) * ...
	(9566.67+-1.59*t5+0.0422*t5*t5) + 1.0 * x5_1 * (-442.13 + v4 ) ...
	*(0.0+15.97*t5+0.0422*t5*t5) + ...
	(v4*x4_1*-1.0*8.547008547008547e-4*exp(18.5751+(-3632.649/(t4+239.2)))) * ...
	(9566.67+-1.59*t4+0.0422*t4*t4) + -1.0 * x6_1 * (-442.13 + v5 ) ...
	*(0.0+15.97*t6+0.0422*t6*t6) + ...
	(v5*x5_2*1.0*8.620689655172415e-4*exp(18.3443+(-3841.2203/(t5+228.0)))) * ...
	(10834.67+8.74*t5+0.0*t5*t5) + 1.0 * x5_2 * (-442.13 + v4 ) ...
	*(0.0+18.1*t5+0.0*t5*t5) + ...
	(v4*x4_2*-1.0*8.547008547008547e-4*exp(18.3443+(-3841.2203/(t4+228.0)))) * ...
	(10834.67+8.74*t4+0.0*t4*t4) + -1.0 * x6_2 * (-442.13 + v5 ) ...
	*(0.0+18.1*t6+0.0*t6*t6))/100000.0))^2 + ...
	((((v6*x6_1*1.0*8.695652173913045e-4*exp(18.5751+(-3632.649/(t6+239.2)))) * ...
	(9566.67+-1.59*t6+0.0422*t6*t6) + 1.0 * x6_1 * (-442.13 + v5 ) ...
	*(0.0+15.97*t6+0.0422*t6*t6) + ...
	(v5*x5_1*-1.0*8.620689655172415e-4*exp(18.5751+(-3632.649/(t5+239.2)))) * ...
	(9566.67+-1.59*t5+0.0422*t5*t5) + -1.0 * x7_1 * (-442.13 + v6 ) ...
	*(0.0+15.97*t7+0.0422*t7*t7) + ...
	(v6*x6_2*1.0*8.695652173913045e-4*exp(18.3443+(-3841.2203/(t6+228.0)))) * ...
	(10834.67+8.74*t6+0.0*t6*t6) + 1.0 * x6_2 * (-442.13 + v5 ) ...
	*(0.0+18.1*t6+0.0*t6*t6) + ...
	(v5*x5_2*-1.0*8.620689655172415e-4*exp(18.3443+(-3841.2203/(t5+228.0)))) * ...
	(10834.67+8.74*t5+0.0*t5*t5) + -1.0 * x7_2 * (-442.13 + v6 ) ...
	*(0.0+18.1*t7+0.0*t7*t7))/100000.0))^2;

        lb = -inf;
        ub = inf; 

    case 145
        %nasty  12949826
        n = 2;
        x0 = [1.0e-30 1.0]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = 0.5*(1.0e10*x1)*(1.0e10*x1) + ...
            0.5*(x2)*(x2);

        lb = -inf;
        ub = inf;   

    case 146
        %osborneb   12949840
        n = 11;
        m = 65;
        x0 = [1.3 0.65 0.65 0.7 0.6 ...
              3 5 7 2 4.5 5.5]';

        % Define the objective function
        x = optimvar('x',n,1);

        y = [1.366 1.191 1.112 1.013 0.991 ...
             0.885 0.831 0.847 0.786 0.725 ...
             0.746 0.679 0.608 0.655 0.616 ...
             0.606 0.602 0.626 0.651 0.724 ...
             0.649 0.649 0.694 0.644 0.624 ...
             0.661 0.612 0.558 0.533 0.495 ...
             0.500 0.423 0.395 0.375 0.372 ...
             0.391 0.396 0.405 0.428 0.429 ...
             0.523 0.562 0.607 0.653 0.672 ...
             0.708 0.633 0.668 0.645 0.632 ...
             0.591 0.559 0.597 0.625 0.739 ...
             0.710 0.729 0.720 0.636 0.581 ...
             0.428 0.292 0.162 0.098 0.054]';

        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        x5 = x(5);
        x6 = x(6);
        x7 = x(7);
        x8 = x(8);
        x9 = x(9);
        x10 = x(10);
        x11 = x(11);

        var = (1:m)';
        t = (var-1)/10;

        objfun = sum( ...
            ( y-x1*exp(-t*x5) - ...
            x2*exp(-((t-x9).^2)*x6) - ...
            x3*exp(-((t-x10).^2)*x7) - ...
        	x4*exp(-((t-x11).^2)*x8) ).^2 ...
            );

        lb = -inf;
        ub = inf;   
   
    case 147
        %yfitu  12951766 
        n = 3;
        x0 = [0.6 -0.6 20]';

        % Define the objective function
        x = optimvar('x',n,1);

        alpha = x(1);
        beta = x(2);
        dist = x(3);

        objfun = (dist * (tan(alpha*(1.0-(0.0/16.0))+beta*(0.0/16.0))) - 21.158931)*(dist * ...
	(tan(alpha*(1.0-(0.0/16.0))+beta*(0.0/16.0))) - 21.158931) + (dist * ...
	(tan(alpha*(1.0-(1.0/16.0))+beta*(1.0/16.0))) - 17.591719)*(dist * ...
	(tan(alpha*(1.0-(1.0/16.0))+beta*(1.0/16.0))) - 17.591719) + (dist * ...
	(tan(alpha*(1.0-(2.0/16.0))+beta*(2.0/16.0))) - 14.046854)*(dist * ...
	(tan(alpha*(1.0-(2.0/16.0))+beta*(2.0/16.0))) - 14.046854) + (dist * ...
	(tan(alpha*(1.0-(3.0/16.0))+beta*(3.0/16.0))) - 10.519732)*(dist * ...
	(tan(alpha*(1.0-(3.0/16.0))+beta*(3.0/16.0))) - 10.519732) + (dist * ...
	(tan(alpha*(1.0-(4.0/16.0))+beta*(4.0/16.0))) - 7.0058392)*(dist * ...
	(tan(alpha*(1.0-(4.0/16.0))+beta*(4.0/16.0))) - 7.0058392) + (dist * ...
	(tan(alpha*(1.0-(5.0/16.0))+beta*(5.0/16.0))) - 3.5007293)*(dist * ...
	(tan(alpha*(1.0-(5.0/16.0))+beta*(5.0/16.0))) - 3.5007293) + (dist * ...
	(tan(alpha*(1.0-(6.0/16.0))+beta*(6.0/16.0))))*(dist * ...
	(tan(alpha*(1.0-(6.0/16.0))+beta*(6.0/16.0)))) + (dist * ...
	(tan(alpha*(1.0-(7.0/16.0))+beta*(7.0/16.0))) + 3.5007293)*(dist * ...
	(tan(alpha*(1.0-(7.0/16.0))+beta*(7.0/16.0))) + 3.5007293) + (dist * ...
	(tan(alpha*(1.0-(8.0/16.0))+beta*(8.0/16.0))) + 7.0058392)*(dist * ...
	(tan(alpha*(1.0-(8.0/16.0))+beta*(8.0/16.0))) + 7.0058392) + (dist * ...
	(tan(alpha*(1.0-(9.0/16.0))+beta*(9.0/16.0))) + 10.519732)*(dist * ...
	(tan(alpha*(1.0-(9.0/16.0))+beta*(9.0/16.0))) + 10.519732) + (dist * ...
	(tan(alpha*(1.0-(10.0/16.0))+beta*(10.0/16.0))) + 14.046854)*(dist * ...
	(tan(alpha*(1.0-(10.0/16.0))+beta*(10.0/16.0))) + 14.046854) + (dist * ...
	(tan(alpha*(1.0-(11.0/16.0))+beta*(11.0/16.0))) + 17.591719)*(dist * ...
	(tan(alpha*(1.0-(11.0/16.0))+beta*(11.0/16.0))) + 17.591719) + (dist * ...
	(tan(alpha*(1.0-(12.0/16.0))+beta*(12.0/16.0))) + 21.158931)*(dist * ...
	(tan(alpha*(1.0-(12.0/16.0))+beta*(12.0/16.0))) + 21.158931) + (dist * ...
	(tan(alpha*(1.0-(13.0/16.0))+beta*(13.0/16.0))) + 24.753206)*(dist * ...
	(tan(alpha*(1.0-(13.0/16.0))+beta*(13.0/16.0))) + 24.753206) + (dist * ...
	(tan(alpha*(1.0-(14.0/16.0))+beta*(14.0/16.0))) + 28.379405)*(dist * ...
	(tan(alpha*(1.0-(14.0/16.0))+beta*(14.0/16.0))) + 28.379405) + (dist * ...
	(tan(alpha*(1.0-(15.0/16.0))+beta*(15.0/16.0))) + 32.042552)*(dist * ...
	(tan(alpha*(1.0-(15.0/16.0))+beta*(15.0/16.0))) + 32.042552) + (dist * ...
	(tan(alpha*(1.0-(16.0/16.0))+beta*(16.0/16.0))) + 35.747869)*(dist * ...
	(tan(alpha*(1.0-(16.0/16.0))+beta*(16.0/16.0))) + 35.747869); 

        lb = -inf;
        ub = inf;

    case 148
        %denschna   12951787 
        n = 2;
        x0 = [1 1]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = x1^4 + (x1+x2)^2 + ...
            (-1.0+exp(x2))^2;

        lb = -inf;
        ub = inf; 

    case 149
        %denschnb   12951819
        n = 2;
        x0 = [1 1]';

        % Define the objective function
        x = optimvar('x',n,1);

        x1 = x(1);
        x2 = x(2);

        objfun = (x1-2.0)^2 + ((x1-2.0)*x2)^2 + ...
            (x2+1.0)^2;

        lb = -inf;
        ub = inf; 

    case 150
        %dixmaank   其余的仅参数不同 12952012
        m = 1000;
        n = 3*m;
        x0 = 2*ones(n,1);

        % Define the objective function
        x = optimvar('x',n,1);

        alpha = 1.0;
        beta = 0.125;
        gamma = 0.125;
        delta  = 0.125;
        K = [2 0 0 2]';

        var1 = (1:n)';
        var2 = (1:n-1)';
        var3 = (1:2*m)';
        var4 = (1:m)';

        xTerm1 = x(1:n-1);  %1,n-1
        xTerm2 = x(2:n);    %2,n
        xTerm3 = x(1:2*m);
        xTerm4 = x((1+m):(2*m+m));
        xTerm5 = x(1:m);
        xTerm6 = x((1+2*m):(m+2*m));

        objfun = 1 + ...
                sum( alpha* ((x.^2).*((var1/n).^K(1))) ) + ...
                sum( beta* ((xTerm1.^2).*((xTerm2+(xTerm2.^2)).^2) .* ((var2/n).^K(2))) )+ ...
                sum( gamma* ((xTerm3.^2) .* (xTerm4.^4) .* ((var3/n).^K(3))) ) + ...
                sum( delta* (xTerm5.*xTerm6.*((var4/n).^K(4))) );

        lb = -inf;
        ub = inf;         

%     case 151          %自动微分不支持变量作为指数
%         %pfit3
%         n = 3;
%         x0 = [1 0 1]';
% 
%         % Define the objective function
%         x = optimvar('x',n,1);
% 
%         a = x(1);
%         r = x(2);
%         h = x(3);
% 
%         cf = -56-(8/9);
%         cg = -126-(2/9);
%         ch = -143-(407/999);
%         
%         part1 = (-0.5*(a*(a+1)*r*h^2)+a*r*h-r*(1-(1+h)^(-a))-cf)^2;
%         part2 = (-a*(a+1)*r*h^2+a*r*h*(1-(1+h)^(-(a+1)))-cg)^2;
%         part3 = (-a*(a+1)*r*h^2*(1-(1+h)^(-(a+2)))-ch)^2;
% 
%         objfun = part1 + part2 + part3;
% 
%         lb = -inf;
%         ub = inf; 
end

% We call the subroutine function CNMGE.m to get the global optimal 
% solution of the problem, the number of local optimal solutions obtained 
% by CNMGE and the computational time.

% The analytical gradient of the objective function is given. Then, call
% this subroute as follows. 
% [x_opt,f_opt,Num_stp, CPUT_time] = CNMGE(func,x0, ub, lb);

% It needs to compute the numerical gradient of the objective function.
% Then, call this subroute as follows.

% In order to generate the correct objective function and gradient function
% by the automatic differentiation fuction prob2strct.m, we add a simple
% constraint as follows. 
% unitdisk = x(1)^2 <= 1; 

addpath(genpath(pwd));
CNMGE_path = which('CNMGE');
CNMGE_path(end-6:end) = [];

unitdisk = x(1)^2 + x(2)^3 <= 1; 
prob = optimproblem("Objective",objfun,"Constraints",unitdisk);
problem = prob2struct(prob,"ObjectiveDerivative","auto-reverse",'FileLocation',CNMGE_path);

%% test_CNMGE
[f_opt, x_opt,CPU_time,Num_stp] = CNMGE(problem.objective, x0, ub, lb);

% Output the global minumum value computed by CNMGE.

fprintf('%%%%%%%%%% %d\n',test);

fprintf('The global minimum value computed by CNMGE is %12.8f\n',f_opt);

%Output the computational time of CNMGE.

fprintf('The computational time of CNMGE is %8.4f\n',CPU_time);

