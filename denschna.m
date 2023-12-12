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
% #   Source: an example problem (p. 206) in
% #   J.E. Dennis and R.B. Schnabel,
% #   "Numerical Methods for Unconstrained Optimization and Nonlinear
% #   Equations",
% #   Prentice-Hall, Englewood Cliffs, 1983.
%
% #   SIF input: Ph. Toint, Nov 1990.
%
% #   classification OUR2-AN-2-0
%
function func = denschna(x)
%denschna   12951787
n = 2;

x1 = x(1);
x2 = x(2);

func =  x1^4 + (x1+x2)^2 + ...
    (-1.0+exp(x2))^2;

end

