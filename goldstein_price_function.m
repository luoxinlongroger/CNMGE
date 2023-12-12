% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = goldstein_price_function(x);
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
function [func,grad_func] = goldstein_price_function(x)
x1 = x(1);
x2 = x(2);
func = (1+(x1+x2+1)^2*(19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2))*...
 (30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2));

% Compute the gradient of the objective function at x. 
if nargout > 1
     x1 = x(1);
     x2 = x(2);
     n = length(x);
     grad_func = zeros(n,1);
     z1 = (1+(x1+x2+1)^2*(19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2));
     z2 = (30+(2*x1-3*x2)^2*(18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2));
     p2 = (19-14*x1+3*x1^2-14*x2+6*x1*x2+3*x2^2);
     p4 = (18-32*x1+12*x1^2+48*x2-36*x1*x2+27*x2^2);
     grad_func(1) = (2*(x1+x2+1)*p2)*z2 +...
         ((x1+x2+1)^2*(-14+6*x1+6*x2))*z2 + z1*(4*(2*x1-3*x2)*p4) +...
         z1*((2*x1-3*x2)^2*(-32+24*x1-36*x2));
     grad_func(2) = (2*(x1+x2+1)*p2)*z2 +...
         ((x1+x2+1)^2*(-14+6*x1+6*x2))*z2 + z1*(-6*(2*x1-3*x2)*p4) +...
         z1*((2*x1-3*x2)^2*(48-36*x1+54*x2));
end

end

