% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = raydan_2_function(x);
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
function [func,grad_func] = extended_quadratic_penalty_qp2_function(x)
n = length(x);
sumx2 = x'*x - 100;
func = 0;
for i = 1 : n-1
    func = func + (x(i)^2-sin(x(i)))^2;
end
func = func + sumx2^2;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    sumx2 = x'*x - 100;
    for i = 1 : n-1
        grad_func(i) = 2*(x(i)^2-sin(x(i)))*(2*x(i)-cos(x(i))) + 4*x(i)*sumx2;
    end
    grad_func(n) = 4*x(n)*sumx2;
end 

end

