% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = himmelblau_function(x);
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
function [func,grad_func] = himmelblau_function(x)
x1 = x(1);
x2 = x(2);
func = (x1^2 + x2 -11)^2 + (x1 + x2^2 -7)^2;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    x1 = x(1);
    x2 = x(2);
    grad_func = zeros(n,1);
    grad_func(1) = 4*(x1^2 + x2 -11)*x1 + 2*(x1 + x2^2 -7);
    grad_func(2) = 2*(x1^2 + x2 -11) + 4*(x1 + x2^2 -7)*x2;
end 

end

