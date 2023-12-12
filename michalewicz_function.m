% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = michalewicz_function(x);
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
% grad_func: It denotes the gradient of the objective function at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function [func,grad_func] = michalewicz_function(x)
n = length(x);
func = 0;
for i = 1 : n
    func = func - sin(x(i))*(sin(i*x(i)^2/pi)^20);
end

% Compute the gradient of the objective function at x.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n
        grad_func(i) = -cos(x(i))*(sin(i*x(i)^2/pi)^20)-...
            40*i*x(i)*sin(x(i))*(sin((i*x(i)^2)/pi)^19)*cos((i*x(i)^2)/pi)/pi;
    end
end

end

