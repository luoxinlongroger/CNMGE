%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PERM FUNCTION d, beta
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
% xx = [x1, x2]
% b  = constant (optional), with default value 0.5
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = perm_function(x);
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
% func: It denotes the objective function. 
% 
% grad_func: It denotes the gradient of the objective function at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
%
function [func,grad_func] = perm_function(x)
b = 10;
n = length(x);
outer = 0;

for ii = 1:n
	inner = 0;
	for jj = 1:n
		xj = x(jj);
        inner = inner + (jj^ii+b)*((xj/jj)^ii-1);
	end
	outer = outer + inner^2;
end

func = outer;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    beta = 10;
    temp1 = 0;
    temp2 = 0;
    for k = 1 :  n
        for i = 1 : n
            for j = 1 : n
                temp1 = temp1 + (j^i+beta)*((x(j)/j)^i-1);
            end
            temp2 = temp2 + i*(x(k)/k)^(i-1)*(k^i+beta)*temp1/k;
            temp1 = 0;
        end
        grad_func(k) = 2*temp2;
        temp2 = 0;
    end
end

end

