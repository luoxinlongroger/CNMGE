% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = hosaki_function(x);
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
% [1] Adorio, E. P., Diliman, U. P.: MVF-Multivariate test functions 
% library in C for unconstrained global optimization, available at
% http://www.geocities.ws/eadorio/mvf.pdf, 2005.
%
function [func,grad_func] = hosaki_function(x)
func = (1-8*x(1)+7*x(1)^2-(7*x(1)^3)/3+(x(1)^4)/4)*(x(2)^2)*exp(-x(2));

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = (-8+14*x(1)-7*x(1)^2+x(1)^3)*(x(2)^2)*exp(-x(2));
    grad_func(2) = (1-8*x(1)+7*x(1)^2-(7*x(1)^3)/3+(x(1)^4)/4)*...
        (2*x(2)*exp(-x(2))-(x(2)^2)*exp(-x(2)));
end 

end

