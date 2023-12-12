%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% COLVILLE FUNCTION
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
% x = [x1, x2, x3, x4]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = colville_function(x)
x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

term1 = 100*(x1^2-x2)^2;
term2 = (x1-1)^2;
term3 = (x3-1)^2;
term4 = 90*(x3^2-x4)^2;
term5 = 10.1*((x2-1)^2 + (x4-1)^2);
term6 = 19.8*(x2-1)*(x4-1);

func = term1 + term2 + term3 + term4 + term5 + term6;

% Compute the gradient of the objective function.
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    grad_func(1) = 400*x(1)*(x(1)^2-x(2))+2*(x(1)-1);
    grad_func(2) = -200*(x(1)^2-x(2))+20.2*(x(2)-1)+19.8*(x(4)-1);
    grad_func(3) = 2*(x(3)-1)+360*x(3)*(x(3)^2-x(4));
    grad_func(4) = -180*(x(3)^2-x(4))+20.2*(x(4)-1)+19.8*(x(2)-1);
end

end

