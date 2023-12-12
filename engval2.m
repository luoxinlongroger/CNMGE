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
% #   Source: problem 15 in
% #   A.R. Buckley,
% #   "Test functions for unconstrained minimization",
% #   TR 1989CS-3, Mathematics, statistics and computing centre,
% #   Dalhousie University, Halifax (CDN), 1989.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-3-0
%
function func = engval2(x)
%engval2        12882702
n = 3;

x1 = x(1);
x2 = x(2);
x3 = x(3);

func = (x1^2+x2^2+x3^2-1)^2   ...
    + (x1^2+x2^2+(x3-2)^2-1)^2  ...
    + (x1+x2+x3-1)^2    ...
    + (x1+x2-x3+1)^2    ....
    + (3*x2^2+x1^3+(5*x3-x1+1)^2-36)^2;


end

