% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = booth_function(x);
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
function [func,grad_func] = drop_wave_function(x)
func = -(1+cos(12*sqrt(x(1)^2+x(2)^2)))/(0.5*(x(1)^2+x(2)^2)+2);

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    a = x(1)^2+x(2)^2;
    sqa = sqrt(a);
    sina = sin(12*sqa);
    cosa = cos(12*sqa);
    grad_func(1) = (sina*12*x(1)/sqa)/(0.5*a+2)+(1+cosa)*x(1)/((0.5*a+2)^2);
    grad_func(2) = (sina*12*x(2)/sqa)/(0.5*a+2)+(1+cosa)*x(2)/((0.5*a+2)^2);

end 

end

