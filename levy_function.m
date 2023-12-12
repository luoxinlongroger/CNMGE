%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LEVY FUNCTION
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
% INPUT:
%
% x = [x1, x2, ..., xd]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = levy_function(x);
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
% func: It denotes the objective function at x. 
%
% grad_func: It denotes the gradient of the objective function at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
function [func,grad_func] = levy_function(x)
n = length(x);

w = 1 + (x-1)/4;
term1 = (sin(pi*w(1)))^2;
term3 = (w(n)-1)^2 * (1+(sin(2*pi*w(n)))^2);

term2 = 0;
for ii = 1:(n-1)
	wi = w(ii);
    new = (wi-1)^2 * (1+10*(sin(pi*wi+1))^2);
	term2 = term2 + new;
end

func = term1 + term2 + term3;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    w = 1+(x-1)/4;
    grad_func(1) = 2*pi*sin(pi*w(1))*cos(pi*w(1))/4+...
        (0.5)*(w(1)-1)*(1+10*sin(pi*w(1)+1)^2)+...
        (w(1)-1)^2*(5*pi*sin(pi*w(1)+1)*cos(pi*w(1)+1));
    for i = 2 : n-1
        grad_func(i) = (0.5)*(w(i)-1)*(1+10*sin(pi*w(i)+1)^2)+...
            (w(i)-1)^2*(5*pi*sin(pi*w(i)+1)*cos(pi*w(i)+1));
    end
    grad_func(n) = 0.5*(w(n)-1)*(1+(sin(2*pi*w(n)))^2)+...
        (w(n)-1)^2*(pi*sin(2*pi*w(n))*cos(2*pi*w(n)));
end

end

