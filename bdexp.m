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
% #   Source: Problem 56 in
% #   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
% #   "Performance of a multifrontal scheme for partially separable
% #   optimization",
% #   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.
% 
% #   SIF input: Ph. Toint, Dec 1989.
% 
% #   classification OBR2-AY-V-0
%
function func = bdexp(x)
%bdexp
n = 1000;
ngs = n-2;

term1 = x(1:ngs);
term2 = x(2:ngs+1);
term3 = x(3:ngs+2);

func = sum( (term1 + term2) ...
    .* exp((term1+term2).*(-term3)) );


end

