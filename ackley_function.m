% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = ackley_function(x);
%
% This is a test problem subroutine used to test the performance of CNMGE.
% The test problem comes from [1].
%
% Input:
% x: It denotes the current x.
%
% Output: 
% func: It denotes the function value at x.
%
% grad_func: It denotes the gradient of the objective function  at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function [func,grad_func] = ackley_function(x)
d = 1000;
c = 2*pi;
b = 0.2;
a = 20;
sum1 = x'*x;
cosy = cos(c*x);
sum2 = sum(cosy);
term1 = -a * exp(-b*sqrt(sum1/d));
term2 = -exp(sum2/d);
func = term1 + term2 + a + exp(1);
% Output the gradient of the objective function.
if nargout > 1
    grad_func = -a*exp(-b*sqrt(sum1/d))*(-b*x/(sqrt(d*sum1)))...
        -exp(sum2/d)*(-c/d*sin(c*x));
end

end

