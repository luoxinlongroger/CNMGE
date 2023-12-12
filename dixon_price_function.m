%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DIXON-PRICE FUNCTION
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = dixon_price_function(x)
x1 = x(1);
n = length(x);
term1 = (x1-1)^2;
sumfun = 0;
for ii = 2:n
	xi = x(ii);
	xold = x(ii-1);
	new = ii * (2*xi^2 - xold)^2;
	sumfun = sumfun + new;
end
func = term1 + sumfun;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = 2*(x(1)-1)-4*(2*x(2)^2-x(1));
    for i = 2 : n-1
        grad_func(i) = 8*i*x(i)*(2*x(i)^2-x(i-1))-2*(i+1)*(2*x(i+1)^2-x(i));
    end
    grad_func(n) = 8*n*x(n)*(2*x(n)^2-x(n-1));
end 

end

