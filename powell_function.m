%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% POWELL FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUTS:
%
% x = [x1, x2, ..., xn]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func, grad_func] = powell_function(x);
%
% This is a subroutine that records the gradient of the test problem and
% used to test the performance of CNMGE.
% The test problem comes from [1], and the gradient of the test problem is
% calculated by Hang Xiao and Xinlong Luo.
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
function [func,grad_func] = powell_function(x)
n = length(x);
sumfun = 0;
for ii = 1:(n/4)
	term1 = (x(4*ii-3) + 10*x(4*ii-2))^2;
	term2 = 5 * (x(4*ii-1) - x(4*ii))^2;
	term3 = (x(4*ii-2) - 2*x(4*ii-1))^4;
	term4 = 10 * (x(4*ii-3) - x(4*ii))^4;
	sumfun = sumfun + term1 + term2 + term3 + term4;
end
func = sumfun;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n/4
        grad_func(4*i-3) = 2*(x(4*i-3)+10*x(4*i-2))+40*(x(4*i-3)-x(4*i))^3;
        grad_func(4*i-2) = 20*(x(4*i-3)+10*x(4*i-2))+4*(x(4*i-2)-2*x(4*i-1))^3;
        grad_func(4*i-1) = 10*(x(4*i-1)-x(4*i))-8*(x(4*i-2)-2*x(4*i-1))^3;
        grad_func(4*i) = -10*(x(4*i-1)-x(4*i))-40*(x(4*i-3)-x(4*i))^3;  
    end
end

end

