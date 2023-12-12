%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EGGHOLDER FUNCTION
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
% x = [x1, x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = eggholder_function(x)
x1 = x(1);
x2 = x(2);

term1 = -(x2+47) * sin(sqrt(abs(x2+x1/2+47)));
term2 = -x1 * sin(sqrt(abs(x1-(x2+47))));

func = term1 + term2;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    x1 = x(1);
    x2 = x(2)+47;
    flag1 = x1/2+x2;
    flag2 = x1-x2;
    sqr1 = sqrt(abs(flag1));
    sqr2 = sqrt(abs(flag2));
    sin1 = sin(sqr1);
    sin2 = sin(sqr2);
    cos1 = cos(sqr1);
    cos2 = cos(sqr2);
    grad_func = zeros(n,1);
    if(flag1 >= 0 && flag2 >= 0)
        grad_func(1) = -x2*cos1/(4*sqr1)-sin2-x1*cos2/(2*sqr2);
        grad_func(2) = -sin1-x2*cos1/(2*sqr1)+x1*cos2/(2*sqr2);
    elseif(flag1 >= 0 && flag2 < 0)
        grad_func(1) = -x2*cos1/(4*sqr1)-sin2+x1*cos2/(2*sqr2);
        grad_func(2) = -sin1-x2*cos1/(2*sqr1)-x1*cos2/(2*sqr2);
    elseif(flag1 < 0 && flag2 >=0 )
        grad_func(1) = x2*cos1/(4*sqr1)-sin2-x1*cos2/(2*sqr2);
        grad_func(2) = -sin1+x2*cos1/(2*sqr1)+x1*cos2/(2*sqr2);
    else
        grad_func(1) = x2*cos1/(4*sqr1)-sin2+x1*cos2/(2*sqr2);
        grad_func(2) = -sin1+x2*cos1/(2*sqr1)-x1*cos2/(2*sqr2);
    end
end 

end

