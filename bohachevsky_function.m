% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = bohachevsky_function(x);
%
% This is a test problem subroutine used to test the performance of CNMGE.
% The test problem comes from [1].
%
% Input:
% x: It denotes the current x.
%
% Output: 
%
% func: It denotes the function value at x.
%
% grad_func: It denotes the gradient of the objective function  at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function [func,grad_func] = bohachevsky_function(x)
func = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))-0.4*cos(4*pi*x(2))+0.7;

% Compute the gradient of the objective function. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = 2*x(1)+0.9*pi*sin(3*pi*x(1));
    grad_func(2) = 4*x(2)+1.6*pi*sin(4*pi*x(2));
end 

end

