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
% #   Source:  problem 8 in
% #   A.R.Conn,N.I.M.Gould and Ph.L.Toint,
% #   "Testing a class of methods for solving minimization
% #   problems with simple bounds on their variables,
% #   Mathematics of Computation 50, pp 399-430, 1988.
%
% #   SIF input: Nick Gould and Ph. Toint, Dec 1995.
%
% #   classification SUR2-AN-V-0
%
function func = chainwoo(x)
%chainwoo   12872931
ns = 499;
n = 2*ns + 2;

term1 = x(1:2:n-2); %x[2*i-1]
term2 = x(2:2:n-2); %x[2*i]
term3 = x(3:2:n);   %x[2*i+1]
term4 = x(4:2:n);   %x[2*i+2]


func = 1 + ...
    sum ( 100*((term2-term1.^2).^2) + ...
    ((1.0-term1).^2) + ...
    90*((term4-term3.^2).^2) + ...
    ((1.0-term3).^2) + ...
    10*((term2+term4-2.0).^2) + ...
    (((term2-term4).^2)./10) );


end

