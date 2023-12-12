%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HARTMANN 3-DIMENSIONAL FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%
%          Derek Bingham, Simon Fraser University
%
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
%
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
%
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
% xx = [x1, x2, x3]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = hartmann_3_d_dimensional_function(x)
alpha = [1.0, 1.2, 3.0, 3.2]';

A = [3.0, 10, 30;
     0.1, 10, 35;
     3.0, 10, 30;
     0.1, 10, 35];

P = 10^(-4) * [3689, 1170, 2673;
               4699, 4387, 7470;
               1091, 8732, 5547;
               381, 5743, 8828];
outer = 0;
for ii = 1:4
	inner = 0;
    for jj = 1:3
        xj = x(jj);
        Aij = A(ii, jj);
        Pij = P(ii, jj);
        inner = inner + Aij*(xj-Pij)^2;
    end
	new = alpha(ii) * exp(-inner);
	outer = outer + new;
end

func = -outer;

% Compute the gradient of the objective function at x. 
if nargout > 1
    alpha = [1.0, 1.2, 3.0, 3.2]';

    A = [3.0, 10, 30;
         0.1, 10, 35;
         3.0, 10, 30;
         0.1, 10, 35];
    
    P = 10^(-4) * [3689, 1170, 2673;
                   4699, 4387, 7470;
                   1091, 8732, 5547;
                   381.5, 5743, 8828];
    n = 3;
    grad_func = zeros(n,1);
    
    grad_func(1) = 2*A(1,1)*(x(1)-P(1,1))*alpha(1)*...
        exp(-A(1,1)*(x(1)-P(1,1))^2-A(1,2)*(x(2)-P(1,2))^2-...
        A(1,3)*(x(3)-P(1,3))^2)+2*A(2,1)*(x(1)-P(2,1))*alpha(2)*...
        exp(-A(2,1)*(x(1)-P(2,1))^2-A(2,2)*(x(2)-P(2,2))^2-A(2,3)*...
        (x(3)-P(2,3))^2)+2*A(3,1)*(x(1)-P(3,1))*alpha(3)*...
        exp(-A(3,1)*(x(1)-P(3,1))^2-A(3,2)*(x(2)-P(3,2))^2-...
        A(3,3)*(x(3)-P(3,3))^2)+2*A(4,1)*(x(1)-P(4,1))*alpha(4)*...
        exp(-A(4,1)*(x(1)-P(4,1))^2-A(4,2)*(x(2)-P(4,2))^2-...
        A(4,3)*(x(3)-P(4,3))^2);
    grad_func(2) = 2*A(1,2)*(x(2)-P(1,2))*alpha(1)*...
        exp(-A(1,1)*(x(1)-P(1,1))^2-A(1,2)*(x(2)-P(1,2))^2-...
        A(1,3)*(x(3)-P(1,3))^2)+2*A(2,2)*(x(2)-P(2,2))*alpha(2)*...
        exp(-A(2,1)*(x(1)-P(2,1))^2-A(2,2)*(x(2)-P(2,2))^2-A(2,3)*...
        (x(3)-P(2,3))^2)+2*A(3,2)*(x(2)-P(3,2))*alpha(3)*...
        exp(-A(3,1)*(x(1)-P(3,1))^2-A(3,2)*(x(2)-P(3,2))^2-...
        A(3,3)*(x(3)-P(3,3))^2)+2*A(4,2)*(x(2)-P(4,2))*alpha(4)*...
        exp(-A(4,1)*(x(1)-P(4,1))^2-A(4,2)*(x(2)-P(4,2))^2-...
        A(4,3)*(x(3)-P(4,3))^2);
    grad_func(3) = 2*A(1,3)*(x(3)-P(1,3))*alpha(1)*...
        exp(-A(1,1)*(x(1)-P(1,1))^2-A(1,2)*(x(2)-P(1,2))^2-...
        A(1,3)*(x(3)-P(1,3))^2)+2*A(2,3)*(x(3)-P(2,3))*alpha(2)*...
        exp(-A(2,1)*(x(1)-P(2,1))^2-A(2,2)*(x(2)-P(2,2))^2-A(2,3)*...
        (x(3)-P(2,3))^2)+2*A(3,3)*(x(3)-P(3,3))*alpha(3)*...
        exp(-A(3,1)*(x(1)-P(3,1))^2-A(3,2)*(x(2)-P(3,2))^2-...
        A(3,3)*(x(3)-P(3,3))^2)+2*A(4,3)*(x(3)-P(4,3))*alpha(4)*...
        exp(-A(4,1)*(x(1)-P(4,1))^2-A(4,2)*(x(2)-P(4,2))^2-...
        A(4,3)*(x(3)-P(4,3))^2);
end

end

