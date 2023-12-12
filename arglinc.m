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
function func = arglinc(x)
%arglinb
n = 10;
m = 20;

% Define the objective function

indexvarN = 1:10;
indexvarN(1) = 0;
indexvarN(10) = 0;
indexvarM = 1:m-1;
indexvarM = indexvarM-1;
indexvarMt = repmat(indexvarM,n,1);
    
func = 2 + sum(((x'.*indexvarN)*indexvarMt-1).^2);


end

