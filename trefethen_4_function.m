% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = trefethen_4_function(x);
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
% [1] Andrei, N.: An unconstrained optimization test functions collection,
% Advanced Modeling and Optimization 10 (2008), 147-161.
%
function [func,grad_func] = trefethen_4_function(x)
func = exp(sin(50*x(1)))+sin(60*exp(x(2)))+sin(70*sin(x(1)))+...
    sin(sin(80*x(2)))-sin(10*(x(1)+x(2)))+0.25*(x(1)^2+x(2)^2);

% Compute the gradient of the objective function at x.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = 50*cos(50*x(1))*exp(sin(50*x(1)))+...
        70*cos(x(1))*cos(70*sin(x(1)))-10*cos(10*(x(1)+x(2)))+x(1)/2;
    grad_func(2) = 60*exp(x(2))*cos(60*exp(x(2)))+...
        80*cos(80*x(2))*cos(sin(80*x(2)))-10*cos(10*(x(1)+x(2)))+x(2)/2;
end 

end

