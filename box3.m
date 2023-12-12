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
% #   Source: Problem 12 in
% #   J.J. More', B.S. Garbow and K.E. Hillstrom,
% #   "Testing Unconstrained Optimization Software",
% #   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
% #   See also Buckley#BOX663
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-3-0
%
function func = box3(x)
%box3   正确  12868149
n = 3;
m = 10;

t = 0.1*(1:m)';
x1 = x(1)*ones(m,1);
x2 = x(2)*ones(m,1);
x3 = x(3)*ones(m,1);

func = sum( (exp(-t.*x1) - exp(-t.*x2) ...
    -x3.*exp(-t) ...
    +x3.*exp(-10*t)).^2 );


end

