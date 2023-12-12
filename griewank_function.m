%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GRIEWANK FUNCTION
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
% x = [x1, x2, ..., xn]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = griewank_function(x)
n = length(x);
sumfun = 0;
prodfun = 1;
for ii = 1:n
	xi = x(ii);
	sumfun = sumfun + xi^2/4000;
	prodfun = prodfun * cos(xi/sqrt(ii));
end
func = sumfun - prodfun + 1;

% Compute the gradient of the objective function at x.
if nargout > 1
    n=length(x);
    a=1;
    sumgrad = zeros(n,1);
    grad_func = zeros(n,1);
    
    for i = 1 : n
        a = a*cos(x(i)/sqrt(i));
    end
    for i = 1 : n
        sumgrad(i) = a*(1/sqrt(i))*sin(x(i)/sqrt(i))/cos(x(i)/sqrt(i));
    end
    for i = 1 : n
        grad_func(i) = x(i)/2000 + sumgrad(i);
    end
end 

end

