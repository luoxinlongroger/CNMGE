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
% This solver finds the global minimum of an unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 
% and quasi-genetic evolutionary algorithms. This program runs all programs 
% successively in desired way.
%
% This is an entrance subroutine when we use this software to find the
% global optimal solution of an unconstrained optimization problem.
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
% func: It denotes the objective function: R^{n}-->R. Since we use automatic
% differentiation to solve its gradient, it includes the objective function 
% and its gradient. If the output of func(x) only has one parameter, it only 
% compute the objective function value at x. Otherwise, it return the objective
% function value and its gradient at x. 
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

        termVer1 = (-2*x(2)+1+(3-2*x(1))*x(1))^(7/3);
        termVer2 = (1-term2 - 2*term3 + (3-2*term1).*term1).^(7/3);
        termVer3 = (-x(n-1)+1+ (3-2*x(n)) *x(n))^(7/3);
        termVer4 = (term4 + term5).^(7/3);
 
        objfun = sqrt(termVer1^2) ...
        + sum( sqrt(termVer2.^2) )...
        + sqrt(termVer3^2)...
        + sum( sqrt(termVer4.^2) );

        lb = -inf;
        ub = inf;

    case 107
        %curly10  12880844
        n = 1000;
        k = 10;
        xi = 1:n;
        x0 = (0.0001*xi/(n+1))';

        % Define the objective function
        x = optimvar('x',n,1);

        X = repmat(x,1,n);
%         Xtr = tril(X);  %下三角 不能用该函数
        e = tril(ones(n));
        Xtr = X.*e;
        XtrNk = [Xtr(:,11:n),zeros(n,10)];

        objfun = sum(sum(Xtr-XtrNk));

        lb = -inf;
        ub = inf;        

        
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

unitdisk = x(1)^2 <= 1; 
prob = optimproblem("Objective",objfun,"Constraints",unitdisk);
problem = prob2struct(prob,"ObjectiveDerivative","auto-reverse",'FileLocation',CNMGE_path);

% test_CNMGE
[f_opt, x_opt,CPU_time,Num_stp] = CNMGE(problem.objective, x0, ub, lb);

% Output the global minumum value computed by CNMGE.

fprintf('The global minimum value computed by CNMGE is %12.8f\n',f_opt);

%Output the computational time of CNMGE.

fprintf('The computational time of CNMGE is %8.4f\n',CPU_time);
