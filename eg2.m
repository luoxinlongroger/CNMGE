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
% #   A.R. Conn, N. Gould and Ph.L. Toint,
% #   "LANCELOT, A Fortran Package for Large-Scale Nonlinear Optimization
% #   (Release A)"
% #   Springer Verlag, 1992.
% 
% #   SIF input: N. Gould and Ph. Toint, June 1994.
% 
% #   classification OUR2-AN-1000-0

%
function func = eg2(x)
%eg2        12882622
n = 1000;

term1 = x(1:n-1);   %xi

func = sum(sin(x(1) + (term1.^2) - 1.0) )...
    + 0.5*sin(x(n)^2);


end

