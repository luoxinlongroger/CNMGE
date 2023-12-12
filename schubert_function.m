% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = schubert_function(x);
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
function [func,grad_func] = schubert_function(x)
n = length(x);
sumfun = 0;
for i = 1 : n
    sumfun = sumfun + sin(2*x(i)+1)+...
        2*sin(3*x(i)+2)+3*sin(4*x(i)+3)+4*sin(5*x(i)+4)+5*sin(6*x(i)+5);
end
func = -sumfun;

% Compute the gradient of the objective function at x.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n
        grad_func(i) = 2*cos(2*x(i)+1)+6*cos(3*x(i)+2)+12*cos(4*x(i)+3)+...
            20*cos(5*x(i)+4)+30*cos(6*x(i)+5);
    end
    grad_func = -grad_func;
end

end

