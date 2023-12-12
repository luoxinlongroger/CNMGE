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
% #   N. Gould, private communication.
% 
% #   SIF input: Nick Gould, June 1990.
% 
% #   classification OUR2-AY-4-0
%
function func = allinitu(x)
%allinitu       12868194
n = 4;

func = x(3)-1 + x(1)^2+ x(2)^2 + (x(3)+x(4))^2 + sin(x(3))^2 + ...
    x(1)^2*x(2)^2 + x(4)-3 + ...
    sin(x(3))^2 + (x(4)-1)^2 + (x(2)^2)^2 +...
    (x(3)^2 + (x(4)+x(1))^2)^2 + (x(1)-4 + ...
    sin(x(4))^2 + x(2)^2*x(3)^2)^2 + sin(x(4))^4;


end

