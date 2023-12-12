function objfleqx = test_function(x)
%FUNCTION 此处显示有关此函数的摘要
%   此处显示详细说明
global test;
switch test
    case 1 %molecular energy problem
        n = 1000;
        k = 4.141720682;
        r = 10.60099896;
        p = 0;
        q = 0;
        for i = 1 : n
            p = p + 1 + cos(3*x(i));
            q = q + ((-1)^i)/sqrt(r-k*cos(x(i)));
        end
        objfleqx = p + q;
    
    case 2 %ackley function
        n = 1000;
        c = 2*pi;
        b = 0.2;
        a = 20;
        sum1 = x'*x;
        cosy = cos(c*x);
        sum2 = sum(cosy);
        term1 = -a * exp(-b*sqrt(sum1/n));
        term2 = -exp(sum2/n);

        objfleqx = term1 + term2 + a + exp(1);
        
    case 3 %levy function
        n = 1000;
        w = 1 + (x-1)/4;
        term1 = (sin(pi*w(1)))^2;
        term3 = (w(n)-1)^2 * (1+(sin(2*pi*w(n)))^2);
        sum1 = 0;
        for ii = 1:(n-1)
            wi = w(ii);
                new = (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
            sum1 = sum1 + new;
        end
        objfleqx = term1 + sum1 + term3;

    case 4 %schwefel function
        n = 1000;
        sum1 = 0;
        for ii = 1:n
            xi = x(ii);
            sum1 = sum1 + xi*sin(sqrt(abs(xi)));
        end
        objfleqx = 418.9829*n - sum1;
        
    case 5 %rastrigin function
        n = 1000;
        sum1 = 0;
        for ii = 1:n
            xi = x(ii);
            sum1 = sum1 + (xi^2 - 10*cos(2*pi*xi));
        end
        objfleqx = 10*n + sum1;
        
    case 6 %styblinski-tang function
        objfleqx = 0;
        n = 1000;
        for i =1 : n
            objfleqx = objfleqx + 0.5*(x(i)^4-16*x(i)^2+5*x(i));
        end
        
    case 7 %trid function
        n = 1000;
        sum1 = (x(1)-1)^2;
        sum2 = 0;

        for ii = 2:n
            xi = x(ii);
            xold = x(ii-1);
            sum1 = sum1 + (xi-1)^2;
            sum2 = sum2 + xi*xold;
        end
        objfleqx = sum1 - sum2;
        
    case 8 %sum squares function
        n = 1000;
        sum1 = 0;
        for ii = 1:n
            xi = x(ii);
            sum1 = sum1 + ii*xi^2;
        end
        objfleqx = sum1;
        
    case 9 %sphere function
        objfleqx = x'*x;
        
    case 10 %rotated hyper-ellipsoid function
        n = 1000;
        outer = 0;
        for ii = 1:n
            inner = 0;
            for jj = 1:ii
                xj = x(jj);
                inner = inner + xj^2;
            end
            outer = outer + inner;
        end
        objfleqx = outer;
        
    case 11 %zakharov function
        n = 1000;
        sum1 = 0;
        sum2 = 0;
        for ii = 1:n
            xi = x(ii);
            sum1 = sum1 + xi^2;
            sum2 = sum2 + 0.5*ii*xi;
        end
        objfleqx = sum1 + sum2^2 + sum2^4;
        
    case 12 %dixon-price function
        x1 = x(1);
        n = 1000;
        term1 = (x1-1)^2;
        sum1 = 0;
        for ii = 2:n
            xi = x(ii);
            xold = x(ii-1);
            new = ii * (2*xi^2 - xold)^2;
            sum1 = sum1 + new;
        end
        objfleqx = term1 + sum1;

    case 13 %rosenbrock function
        n = 1000;
        sum1 = 0;
        for ii = 1:(n-1)
            xi = x(ii);
            xnext = x(ii+1);
            new = 100*(xnext-xi^2)^2 + (xi-1)^2;
            sum1 = sum1 + new;
        end
        objfleqx = sum1;
        
    case 14 %powell function
        n = 1000;
        sum1 = 0;
        for ii = 1:(n/4)
            term1 = (x(4*ii-3) + 10*x(4*ii-2))^2;
            term2 = 5 * (x(4*ii-1) - x(4*ii))^2;
            term3 = (x(4*ii-2) - 2*x(4*ii-1))^4;
            term4 = 10 * (x(4*ii-3) - x(4*ii))^4;
            sum1 = sum1 + term1 + term2 + term3 + term4;
        end
        objfleqx = sum1;
        
    case 15 %quartic with noise function
        if(nargin == 1)
            noise = 0.5;
        end
        objfleqx = 0;
        n = 1000;
        for i = 1 : n
            objfleqx = objfleqx + x(i)^4;
        end
        objfleqx = objfleqx + noise;
        
    case 16 %schubert function
        n = 1000;
        sum1 = 0;
        for i = 1 : n
            sum1 = sum1+sin(2*x(i)+1)+2*sin(3*x(i)+2)+3*sin(4*x(i)+3)+4*sin(5*x(i)+4)+5*sin(6*x(i)+5);
        end
        objfleqx = -sum1;
        
     case 17 %raydan 1 function
        n = 1000;
        objfleqx = 0;
        for i = 1 : n
            objfleqx = objfleqx + (i/10)*(exp(x(i))-x(i));
        end
        
    case 18 %rayden 2 function
        n = 1000;
        expx = exp(x);
        sumexpx = sum(expx);
        sumx = sum(x);
        objfleqx = sumexpx + sumx;
        
    case 19 %extended_tridiagonal_1_function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx + (x(2*i-1) + x(2*i) - 3)^2 + (x(2*i-1) - x(2*i)+1)^4;
        end
        
     case 20 %extended_quadratic_penalty_qp_1_function
        n = 1000;
        sum1 = 0;
        sum2 = (x'*x - 0.5)^2;
        for i = 1 : n-1
            sum1 = sum1 + (x(i)^2 - 2)^2;
        end
        objfleqx = sum1 + sum2;
        
    case 21 %extended_quadratic_penalty_qp_2_function
        n = 1000;
        sum1 = 0;
        sum2 = (x'*x - 100)^2;
        for i = 1 : n-1
            sum1 = sum1 + (x(i)^2 - sin(x(i)))^2;
        end
        objfleqx = sum1 + sum2;
        
    case 22 %quadratic qf2 function
        n = 1000;
        sum1 = 0;
        for i = 1 : n
            sum1 = sum1 + i*(x(i)^2 - 1)^2;
        end
        objfleqx = 0.5*sum1 - x(n);
        
    case 23 %extended psc1 function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx + (x(2*m-1)^2 + x(2*m)^2 + x(2*m-1)*x(2*m))^2 + sin(x(2*m-1))^2 + cos(x(2*m))^2;
        end
        
    case 24 %extended bd1 function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx +  (x(2*i-1)^2 + x(2*i)^2 - 2)^2 + (exp(x(2*i-1)-1)-x(2*i))^2;
        end
        
    case 25 %extended cliff function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx +  ((x(2*i-1)-3)/100)^2 - (x(2*i-1) - x(2*i)) + exp(20*(x(2*i-1)-x(2*i)));
        end
        
    case 26 %perturbed_quadratic_diagonal_function
        n = 1000;
        objfleqx = 0;
        sumx = sum(x);
        for i = 1 : n
            objfleqx = objfleqx + i*x(i)^2/100;
        end
        objfleqx = objfleqx + sumx^2;
        
    case 27 %extended_hiebert_function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx + (x(2*i-1)-10)^2 + (x(2*i-1)*x(2*i)-50000)^2;
        end
        
    case 28 %extended_tet_function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx +  exp(x(2*i-1) + 3*x(2*i)-0.1) + exp(x(2*i-1) - 3*x(2*i)-0.1) + exp(-x(2*i-1)-0.1);
        end
        
    case 29 %diagonal_1_function
        n = 1000;
        objfleqx = 0;
        for i = 1 : n
            objfleqx = objfleqx + exp(x(i)) - i*x(i);
        end
        
    case 30 %diagonal_3_function
        n = 1000;
        objfleqx = 0;
        for i = 1 : n
            objfleqx = objfleqx + exp(x(i)) - i*sin(x(i));
        end
        
    case 31 %diagonal_5_function
        n = 1000;
        objfleqx = 0;
        for i = 1 : n
            objfleqx = objfleqx + log(exp(x(i)) + exp(-x(i)));
        end
        
    case 32 %extended_maratos_function
        n = 1000;
        m = n/2;
        objfleqx = 0;
        for i = 1 : m
            objfleqx = objfleqx + x(2*i-1) + 100*(x(2*i-1)^2 + x(2*i)^2 - 1)^2;
        end
        
    case 33 %eg2_function
        n = 1000;
        objfleqx = 0;
        for i = 1 : n-1
            objfleqx = objfleqx + sin(x(1) + x(i)^2 -1);
        end
        objfleqx = objfleqx + 0.5*sin(x(n)^2);
        
    case 34 %sinquad_function
        n = 1000;
        objfleqx = (x(1)-1)^4 + (x(n)^2-x(1)^2)^2;
        for i = 2 : n-1
            objfleqx = objfleqx + (sin(x(i)-x(n))-x(1)^2 + x(i)^2)^2;
        end
    case 35 %griewank function
        n = 10;
        sum1 = 0;
        prod = 1;
        for ii = 1:n
            xi = x(ii);
            sum1 = sum1 + xi^2/4000;
            prod = prod * cos(xi/sqrt(ii));
        end
        objfleqx = sum1 - prod + 1;
        
    case 36 %levy function n.13
        x1 = x(1);
        x2 = x(2);
        term1 = (sin(3*pi*x1))^2;
        term2 = (x1-1)^2 * (1+(sin(3*pi*x2))^2);
        term3 = (x2-1)^2 * (1+(sin(2*pi*x2))^2);
        objfleqx = term1 + term2 + term3;
        
    case 37 %hosaki function
        objfleqx = (1-8*x(1)+7*x(1)^2-(7*x(1)^3)/3+(x(1)^4)/4)*(x(2)^2)*exp(-x(2));
        
    case 38 %beale function
        x1 = x(1);
        x2 = x(2);
        term1 = (1.5 - x1 + x1*x2)^2;
        term2 = (2.25 - x1 + x1*x2^2)^2;
        term3 = (2.625 - x1 + x1*x2^3)^2;
        objfleqx = term1 + term2 + term3;
        
    case 39 %easom function
        objfleqx = -cos(x(1))*cos(x(2))*exp(-(x(1)-pi)^2-(x(2)-pi)^2);
        
    case 40 %price function
        objfleqx = (2*x(1)^3*x(2)-x(2)^3)^2+(6*x(1)-x(2)^2+x(2))^2;
        
    case 41 %branin function
        x1 = x(1);
        x2 = x(2);
        t = 1 / (8*pi);
        s = 10;
        r = 6;
        c = 5/pi;
        b = 5.1 / (4*pi^2);
        a = 1;
        term1 = a * (x2 - b*x1^2 + c*x1 - r)^2;
        term2 = s*(1-t)*cos(x1);
        objfleqx = term1 + term2 + s;
        
    case 42 %trecanni function
        objfleqx = x(1)^4+4*x(1)^3+4*x(1)^2+x(2)^2;
        
    case 43 %booth function
        objfleqx = (x(1)+2*x(2)-7)^2+(2*x(1)+x(2)-5)^2;
        
    case 44 %matyas function
        objfleqx = 0.26*(x(1)^2+x(2)^2)-0.48*x(1)*x(2);
        
    case 45 %mccormick function
        objfleqx = sin(x(1)+x(2))+(x(1)-x(2))^2-1.5*x(1)+2.5*x(2)+1;
        
    case 46 % power sum function
        n = 4;
        b = [8, 18, 44, 114];
        outer = 0;
        for ii = 1:n
            inner = 0;
            for jj = 1:n
                xj = x(jj);
                inner = inner + xj^ii;
            end
            outer = outer + (inner-b(ii))^2;
        end
        objfleqx = outer;
        
    case 47 %colville function
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        term1 = 100 * (x1^2-x2)^2;
        term2 = (x1-1)^2;
        term3 = (x3-1)^2;
        term4 = 90 * (x3^2-x4)^2;
        term5 = 10.1 * ((x2-1)^2 + (x4-1)^2);
        term6 = 19.8*(x2-1)*(x4-1);
        objfleqx = term1 + term2 + term3 + term4 + term5 + term6;
        
    case 48 %schaffer function n.2
        x1 = x(1);
        x2 = x(2);
        temp1 = (cos(sin(abs(x1^2-x2^2))))^2 - 0.5;
        temp2 = (1 + 0.001*(x1^2+x2^2))^2;
        objfleqx = 0.5 + temp1/temp2;
        
    case 49 %bohachevsky function
        objfleqx = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;
        
    case 50 %three-hump camel function
        objfleqx = 2*x(1)^2-1.05*x(1)^4+x(1)^6/6+x(1)*x(2)+x(2)^2;
        
    case 51 %six-hump camel function
        objfleqx = (4-2.1*x(1)^2+x(1)^4/3)*x(1)^2+x(1)*x(2)+(-4+4*x(2)^2)*x(2)^2;
        
    case 52 %drop-wave function
        objfleqx = -(1+cos(12*sqrt(x(1)^2+x(2)^2)))/(0.5*(x(1)^2+x(2)^2)+2);
        
    case 53 %perm function 0,alpha,beta
        b = 10;
        n = 4;
        outer = 0;
        for ii = 1:n
            inner = 0;
            for jj = 1:n
                xj = x(jj);
                inner = inner + (jj^ii+b)*((xj/jj)^ii-1);
            end
            outer = outer + inner^2;
        end
        objfleqx = outer;
        
    case 54 %hartmann 3-dimensional function
        alpha = [1.0, 1.2, 3.0, 3.2]';
        A = [3.0, 10, 30;
             0.1, 10, 35;
             3.0, 10, 30;
             0.1, 10, 35];
        P = 10^(-4) * [3689, 1170, 2673;
                       4699, 4387, 7470;
                       1091, 8732, 5547;
                       381, 5743, 8828];
        outer = 0;
        for ii = 1:4
            inner = 0;
            for jj = 1:3
                xj = x(jj);
                Aij = A(ii, jj);
                Pij = P(ii, jj);
                inner = inner + Aij*(xj-Pij)^2;
            end
            new = alpha(ii) * exp(-inner);
            outer = outer + new;
        end

        objfleqx = -outer;
        
    case 55 %trefethen4 function
        objfleqx = exp(sin(50*x(1)))+sin(60*exp(x(2)))+sin(70*sin(x(1)))+sin(sin(80*x(2)))-sin(10*(x(1)+x(2)))+0.25*(x(1)^2+x(2)^2);
        
    case 56 %zettl function
        objfleqx = (x(1)^2+x(2)^2-2*x(1))^2+0.25*x(1);
        
    case 57 %exp2 function
        objfleqx = 0;
        for i = 1 : 10
            objfleqx = objfleqx + (exp(-(i-1)*x(1)/10)-5*exp(-(i-1)*(x(2)/10))-exp(-(i-1)/10)+5*exp(-(i-1)))^2;
        end
        
    case 58 %hansen function
        temp1 = 0;
        temp2 = 0;
        for i = 1 : 5
            temp1 = temp1 + i*cos((i+1)*x(2)+i);
            temp2 = temp2 + i*cos((i-1)*x(1)+i);
        end
        objfleqx = temp1*temp2;
        
    case 59 %schaffer function n.4
        x1 = x(1);
        x2 = x(2);
        temp1 = (cos(sin(abs(x1^2-x2^2))))^2 - 0.5;
        temp2 = (1 + 0.001*(x1^2+x2^2))^2;
        objfleqx = 0.5 + temp1/temp2;
        
    case 60 %holder table function
        x1 = x(1);
        x2 = x(2);
        temp1 = sin(x1)*cos(x2);
        temp2 = exp(abs(1 - sqrt(x1^2+x2^2)/pi));
        objfleqx = -abs(temp1*temp2);
        
    case 61 %gramacy&lee(2012) function
        objfleqx = sin(10*pi*x)/(2*x)+(x-1)^4;
        
    case 62 %eggholder function
        x1 = x(1);
        x2 = x(2);

        term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
        term2 = -x1 * sin(sqrt(abs(x1-(x2+47))));

        objfleqx = term1 + term2;
        
    case 63 %michalewitz function
        n = 2;
        objfleqx = 0;
        for i = 1 : n
            objfleqx = objfleqx - sin(x(i))*(sin(i*x(i)^2/pi)^20);
        end
        
    case 64 %box-betts exponential quadratic sum function
        n = 3;
        objfleqx = 0;
        for i = 1 : n
            gx = exp(-0.1*i*x(1))-exp(-0.1*i*x(2))-(exp(-0.1*i)-exp(-i))*x(3);
            objfleqx = objfleqx + gx^2;
        end
        
    case 65 %cross-in-tray function
        x1 = x(1);
        x2 = x(2);
        fact1 = sin(x1)*sin(x2);
        fact2 = exp(abs(100 - sqrt(x1^2+x2^2)/pi));
        objfleqx = -0.0001 * (abs(fact1*fact2)+1)^0.1;
        
    case 66 %himmelblau function
         x1 = x(1);
         x2 = x(2);
         objfleqx = (x1^2 + x2 -11)^2 + (x1 + x2^2 -7)^2;
         
    case 67 %forrester function
         objfleqx = (6*x-2)^2*sin(12*x-4);
         
    case 68 %goldstein-price function
         x1 = x(1);
         x2 = x(2);
         objfleqx = (1+(x1+x2+1)^2*(19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2))*(30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2));

        
end
end

