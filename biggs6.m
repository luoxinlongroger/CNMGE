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
% #   Source: Problem 21 in
% #   A.R. Buckley,
% #   "Test functions for unconstrained minimization",
% #   TR 1989CS-3, Mathematics, statistics and computing centre,
% #   Dalhousie University, Halifax (CDN), 1989.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-6-0
%
function func = biggs6(x)
%biggs6  对了     12868107
n = 6;
m = 13;

x1 = x(1)*ones(m,1);
x2 = x(2)*ones(m,1);
x3 = x(3)*ones(m,1);
x4 = x(4)*ones(m,1);
x5 = x(5)*ones(m,1);
x6 = x(6)*ones(m,1);

func = sum( ( (-exp(-0.1*ones(m,1))) ...
    +(5*exp(-ones(m,1))) ...
    -(3*exp(-0.4*ones(m,1))) ...
    + x3.*exp(-0.1*ones(m,1).*x1) ...
    - x4.*exp(-0.1*ones(m,1).*x2) ...
    + x6.*exp(-0.1*ones(m,1).*x5)).^2 );


end

