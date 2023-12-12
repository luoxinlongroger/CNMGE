% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% [func,grad_func] = gramacy_lee_function(x);
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
function [func,grad_func] = gramacy_lee_function(x)
func = sin(10*pi*x(1))/(2*x(1))+(x(1)-1)^4;

% Compute the gradient of the objective function at x. 
if nargout > 1
    grad_func = 5*pi*cos(10*pi*x(1))/(x(1))-...
        sin(10*pi*x(1))/(2*x(1)^2)+4*(x(1)-1)^3;
end 

end

