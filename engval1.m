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
% #   Source: problem 31 in
% #   Ph.L. Toint,
% #   "Test problems for partially separable optimization and results
% #   for the routine PSPMIN",
% #   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
% 
% #   See also Buckley#172 (p. 52)
% #   SIF input: Ph. Toint and N. Gould, Dec 1989.
% 
% #   classification OUR2-AN-V-0


%
function func = engval1(x)
%engval1     12882629
n = 1000;

term1 = x(1:n-1);   %xi
term2 = x(2:n);     %xi+1

func = sum( ((term1.^2)+(term2.^2)).^2 ) ...
    + sum( (-4*term1+3.0) );


end

