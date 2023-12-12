% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func, grad_func] = raydan_1_function(x);
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
function [func,grad_func] = raydan_1_function(x)
n = length(x);
func = 0;
for i = 1 : n
   func = func + i/10*(exp(x(i))-x(i)); 
end

% Compute the gradient of the objective function at x.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n
       grad_func(i) = i/10*(exp(x(i))-1); 
    end
end 

end

