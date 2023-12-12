% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = quartic_with_noise_function(x,noise);
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
function [func,grad_func] = quartic_with_noise_function(x,noise)
if(nargin == 1)
    noise = 0.5;
end
func = 0;
n = length(x);
for i = 1 : n
    func = func + x(i)^4;
end
func = func + noise;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n
        grad_func(i) = 4*x(i)^3;
    end
end 

end

