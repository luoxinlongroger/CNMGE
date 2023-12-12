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
% #   SIF input: N. Gould, Jan 1996
%
% #   classification OUR2-AN-V-0
%
function func = cosine(x)
%cosine  12878676
n = 1000;

term1 = x(1:n-1);
term2 = x(2:n);

func = sum( cos(-0.5*term2+(term1.^2)) );



end

