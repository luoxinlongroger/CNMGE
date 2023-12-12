% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = styblinski_tang_function(x);
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
function [func,grad_func] = styblinski_tang_function(x)
func = 0;
n = length(x);
for i =1 : n
    func = func + 0.5*(x(i)^4-16*x(i)^2+5*x(i));
end

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1:n
        grad_func(i) = 4*x(i)^3-32*x(i)+5;
    end
end

end

