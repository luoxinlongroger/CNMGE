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
% #   Source: Problem 27 in
% #   J.J. More', B.S. Garbow and K.E. Hillstrom,
% #   "Testing Unconstrained Optimization Software",
% #   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.
% 
% #   See also Buckley#79
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-V-0
%
function func = brownal(x)
%brownal    正确  12868190
n = 10;

term1 = x(1:n-1);
term2 = sum(x);
term3 = (prod(x)-1)^2;

func = sum( (term1 + term2*ones(9,1) - 11).^2 ) + term3;


end

