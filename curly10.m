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

% #   Source: Nick Gould
% 
% #   SIF input: Nick Gould, September 1997.
% 
% #   classification SUR2-AN-V-0
%
function func = curly10(x)
%curly10  12880844
n = 1000;
k = 10;

X = repmat(x,1,n);
%         Xtr = tril(X);  %下三角 不能用该函数
e = tril(ones(n));
Xtr = X.*e;
XtrNk = [Xtr(:,k+1:n),zeros(n,k)];
obj = Xtr-XtrNk;

func = sum(sum(obj));


end

