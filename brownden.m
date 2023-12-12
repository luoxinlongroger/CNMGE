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
% #   Source: Problem 16 in
% #   J.J. More', B.S. Garbow and K.E. Hillstrom,
% #   "Testing Unconstrained Optimization Software",
% #   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
%
% #   See also Buckley#30
% #   SIF input: Ph. Toint, Dec 1989.
%
% #   classification SUR2-AN-4-0
%
function func = brownden(x)
%brownden  12872716
n = 4;
m = 20;

t = 1:m;
tt = (t/5)';
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

func = sum( ((x1+tt*x2-exp(tt)).^2 ...
    + (x3+x4*sin(tt)-cos(tt)).^2).^2 );


end

