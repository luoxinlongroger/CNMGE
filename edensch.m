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
% #   Source:
% #   G. Li,
% #   "The secant/finite difference algorithm for solving sparse
% #   nonlinear systems of equations",
% #   SIAM Journal on Optimization, (to appear), 1990.
% 
% #   SIF input: Ph. Toint, Apr 1990.
% #              minor correction by Ph. Shott, January 1995.
% 
% #   classification OUR2-AN-V-0
%
function func = edensch(x)
%edensch
n = 2000;

term1 = x(1:1:n-1);
term2 = x(2:1:n);

func =  sum( (term1-2).^4 + (term1.*term2- ...
    2*term2).^2 + (term2+1).^2 ) + 16;


end

