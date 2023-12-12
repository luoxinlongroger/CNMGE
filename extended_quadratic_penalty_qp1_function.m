% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function func = extended_quadratic_penalty_qp1_function_r(x);
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
function [func,grad_func] = extended_quadratic_penalty_qp1_function(x)
n = length(x);
sumx2 = x'*x-0.5;
sumfun1 = 0;
for i = 1 : n-1
    sumfun1 = sumfun1 + (x(i)^2-2)^2;
end
func = sumfun1 + sumx2^2;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    sumx2 = x'*x-0.5;
    for i = 1 : n-1
        grad_func(i) = 4*x(i)*(x(i)^2-2) + 4*sumx2*x(i);
    end
    grad_func(n) = 4*x(n)*sumx2;
end 

end

