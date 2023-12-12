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
% #   Source: problem 156 (p. 51) in
% #   A.R. Buckley,
% #   "Test functions for unconstrained minimization",
% #   TR 1989CS-3, Mathematics, statistics and computing centre,
% #   Dalhousie University, Halifax (CDN), 1989.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification QUR2-AN-V-0
%
function func = dixon3dq(x)
%dixon3dq   12868200
n = 10;

term1 = x(2:1:n-1);
term2 = x(3:1:n);

func = (x(1)-1.0)^2 + ...
    sum((term1-term2).^2) + (x(10)-1.0)^2;


end

