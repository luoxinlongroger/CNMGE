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
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function func = bard(x)
%bard
n = 3;
m = 15;

u = 1:m;
v = 16-u;
w = min(u,v);
y = [0.14 0.18 0.22 0.25 0.29 0.32 ...
    0.35 0.39 0.37	0.58 0.73 0.96 ...
    1.34 2.10 4.39];
ut = u';
vt = v';
wt = w';
yt = y';
x1 = x(1) * ones(m,1);
x2 = x(2) * ones(m,1);
x3 = x(3) * ones(m,1);

func = sum((yt - (x1 + ut./(vt.*x2+wt.*x3))).^2);


end

