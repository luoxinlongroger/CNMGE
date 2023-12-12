% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = extended_tridiagonal_1_function(x);
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
function [func,grad_func] = extended_tridiagonal_1_function(x)
n = length(x);
m = n/2;
func = 0;
for i = 1 : m
   func = func + (x(2*i-1)+x(2*i)-1)^2 + (x(2*i-1)-x(2*i)+1)^4;
end

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    m = n/2;
    for i = 1 : m
       grad_func(2*i-1) = 2*(x(2*i-1)+x(2*i)-1) + 4*(x(2*i-1)-x(2*i)+1)^3;
       grad_func(2*i) = 2*(x(2*i-1)+x(2*i)-1) - 4*(x(2*i-1)-x(2*i)+1)^3;
    end
end

end

