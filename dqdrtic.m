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
% #   Source: problem 22 in
% #   Ph. L. Toint,
% #   "Test problems for partially separable optimization and results
% #   for the routine PSPMIN",
% #   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification QUR2-AN-V-0
%
function func = dqdrtic(x)
%dqdrtic    12882615
n = 1000;

term1 = x(1:n-2);   %xi
term2 = x(2:n-1);   %xi+1
term3 = x(3:n);     %xi+2

func = sum((100*(term2.^2)+100*(term3.^2)+(term1.^2)));


end

