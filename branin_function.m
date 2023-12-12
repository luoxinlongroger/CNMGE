%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BRANIN FUNCTION
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
% x = [x1, x2]
% a = constant (optional), with default value 1
% b = constant (optional), with default value 5.1/(4*pi^2)
% c = constant (optional), with default value 5/pi
% r = constant (optional), with default value 6
% s = constant (optional), with default value 10
% t = constant (optional), with default value 1/(8*pi)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func, grad_func] = branin_function(x)
x1 = x(1);
x2 = x(2);
t = 1/(8*pi);
s = 10;
r = 6;
c = 5/pi;
b = 5.1/(4*pi^2);
a = 1;

term1 = a*(x2 - b*x1^2 + c*x1 - r)^2;
term2 = s*(1-t)*cos(x1);

func = term1 + term2 + s;

% Compute the gradient of the objective function 

if nargout > 1
    n = length(x);
    x1 = x(1);
    x2 = x(2);    
    grad_func = zeros(n,1);
    grad_func(1) = 2*a*(x2-b*x1^2+c*x1-r)*(-2*b*x1+c)-s*(1-t)*sin(x1);
    grad_func(2) = 2*a*(x2-b*x1^2+c*x1-r);
end 

end

