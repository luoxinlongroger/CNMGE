% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function function [func,grad_func] = power_sum_function(x);
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
function [func,grad_func] = power_sum_function(x)
n = length(x);
b = [8, 18, 44, 114];
outer = 0;

for ii = 1:n
    inner = 0;
    for jj = 1:n
        xj = x(jj);
        inner = inner + xj^ii;
    end
    outer = outer + (inner-b(ii))^2;
end

func = outer;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    b = [8, 18, 44, 114];
    for i = 1 : n
        grad_func(i) = 2*(x(1)+x(2)+x(3)+x(4)-b(1))+...
            4*x(i)*(x(1)^2+x(2)^2+x(3)^2+x(4)^2-b(2))+...
            6*x(i)^2*(x(1)^3+x(2)^3+x(3)^3+x(4)^3-b(3))+...
            8*x(i)^3*(x(1)^4+x(2)^4+x(3)^4+x(4)^4-b(4));
    end
end

end

