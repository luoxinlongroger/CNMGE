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
% #   Source: Problem 61 in
% #   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
% #   "Performance of a multifrontal scheme for partially separable
% #   optimization",
% #   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification SUR2-AN-V-0
%
function func = bdqrtic(x)
%bdqrtic        12868098
n = 1000;

term1 = x(1:n-4);
term2 = x(2:n-3);
term3 = x(3:n-2);
term4 = x(4:n-1);
func = sum( (-4*term1+3).^2 ) ...
    + sum( (((term1.^2) ...
    + 2*(term2.^2) + 3*(term3.^2) ...
    + 4*(term4.^2) + 5*x(n)^2).^2) );


end

