% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = schaffer_function_n_2(x);
%
% This is a test problem subroutine used to test the performance of CNMGE.
% The test problem comes from [1].
%
% Input:
% x: It denotes the current x.
%
% Output: 
% func: It denotes the objective function value at x.
%
% grad_func: It denotes the gradient of the objective function at x.
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function [func,grad_func] = schaffer_function_n_2(x)
func = 0.5+(sin(sqrt(x(1)^2+x(2)^2))^2-0.5)/((1+0.001*(x(1)^2+x(2)^2))^2);

% Compute the gradient of the objective function at x.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = (2*sin(sqrt(x(1)^2+x(2)^2))*cos(sqrt(x(1)^2+...
        x(2)^2))*x(1))/((1+0.001*(x(1)^2+x(2)^2))^2*sqrt(x(1)^2+...
        x(2)^2))-(0.004*((sin(sqrt(x(1)^2+x(2)^2)))^2-0.5)*x(1))/...
        ((1+0.001*(x(1)^2+x(2)^2))^3);
    grad_func(2)=(2*sin(sqrt(x(1)^2+x(2)^2))*cos(sqrt(x(1)^2+x(2)^2))*x(2))/...
        ((1+0.001*(x(1)^2+x(2)^2))^2*sqrt(x(1)^2+x(2)^2))-...
        (0.004*((sin(sqrt(x(1)^2+x(2)^2)))^2-0.5)*x(2))/...
        ((1+0.001*(x(1)^2+x(2)^2))^3);
end

end

