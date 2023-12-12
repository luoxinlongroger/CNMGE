% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = sinquad_function(x);
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
function [func,grad_func] = sinquad_function(x)
n = length(x);
func = (x(1)-1)^4 + (x(n)^2-x(1)^2)^2;
for i = 2 : n-1
    func = func + (sin(x(i)-x(n))-x(1)^2 + x(i)^2)^2;
end

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
grad_func = zeros(n,1);
sum1 = 0;
sum3 = 0;
for i = 2 : n-1
    sum1 = sum1 + sin(x(i)-x(n))-x(1)^2 + x(i)^2;
    sum3 = sum3 + 2*(sin(x(i)-x(n))-x(1)^2 + x(i)^2)*(-cos(x(i)-x(n)));
end
%少了x(1)
grad_func(1) = 4*(x(1)-1)^3 -4*x(1)*sum1 - 4*(x(n)^2-x(1)^2)*x(1);
for i = 2 : n-1
    grad_func(i) = 2*(sin(x(i)-x(n))-x(1)^2 + x(i)^2)*(cos(x(i)-x(n)) + 2*x(i));
end
grad_func(n) = sum3 + 4*x(n)*(x(n)^2-x(1)^2);

%     n = length(x);
%     grad_func = zeros(n,1);
%     sum1 = 0;
%     sum3 = 0;
%     for i = 2 : n-1
%         sum1 = sum1 + sin(x(i)-x(n))-x(1)^2 + x(i)^2;
%         sum3 = sum3 + 2*(sin(x(i)-x(n))-x(1)^2 + x(i)^2)*(-cos(x(i)-x(n)));
%     end
%     grad_func(1) = 4*(x(1)-1)^3 -4*sum1 - 4*(x(n)^2-x(1)^2)*x(1);
%     for i = 2 : n-1
%         sum2 = 0;
%         for j = 2 : n-1
%             sum2 = sum2 + 2*(sin(x(j)-x(n))-x(1)^2 +...
%                 x(j)^2)*(cos(x(j)-x(n))+2*x(j));
%         end
%         grad_func(i) = sum2;
%     end
%     grad_func(n) = sum3 + 4*x(n)*(x(n)^2-x(1)^2);
end

end

