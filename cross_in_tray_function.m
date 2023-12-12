%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CROSS-IN-TRAY FUNCTION
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [func,grad_func] = cross_in_tray_function(x)
x1 = x(1);
x2 = x(2);

fact1 = sin(x1)*sin(x2);
fact2 = exp(abs(100 - sqrt(x1^2+x2^2)/pi));

func = -0.0001 * (abs(fact1*fact2)+1)^0.1;

% Compute the gradient of the objective function. 
if nargout > 1
        n = length(x);
    grad_func = zeros(n,1);
    x1 = x(1);
    x2 = x(2);
    
    a = 100 - sqrt(x1^2+x2^2)/pi;
    b = sin(x1)*sin(x2);
    d = b*exp(abs(a));
    if(a>0 && b>0)
        grad_func(1) = -0.0001*0.1*(d+1)^(-0.9)*(cos(x1)*sin(x2)*exp(a) +...
            d*(-x1/(pi*sqrt(x1^2+x2^2))));
        grad_func(2) = -0.0001*0.1*(d+1)^(-0.9)*(sin(x1)*cos(x2)*exp(a) +...
            d*(-x2/(pi*sqrt(x1^2+x2^2))));
    end
    if(a>0 && b<0)
        grad_func(1) = -0.0001*0.1*(-d+1)^(-0.9)*(-cos(x1)*sin(x2)*exp(a) +...
            -d*(-x1/(pi*sqrt(x1^2+x2^2))));
        grad_func(2) = -0.0001*0.1*(-d+1)^(-0.9)*(-sin(x1)*cos(x2)*exp(a) +...
            -d*(-x2/(pi*sqrt(x1^2+x2^2))));
    end
    if(a<0 && b>0)
        grad_func(1) = -0.0001*0.1*(d+1)^(-0.9)*(cos(x1)*sin(x2)*exp(-a) +...
            d*(x1/(pi*sqrt(x1^2+x2^2))));
        grad_func(2) = -0.0001*0.1*(d+1)^(-0.9)*(sin(x1)*cos(x2)*exp(-a) +...
            d*(x2/(pi*sqrt(x1^2+x2^2))));
    end
    if(a<0 && b<0)
        grad_func(1) = -0.0001*0.1*(-d+1)^(-0.9)*(-cos(x1)*sin(x2)*exp(-a) +...
            -d*(x1/(pi*sqrt(x1^2+x2^2))));
        grad_func(2) = -0.0001*0.1*(-d+1)^(-0.9)*(-sin(x1)*cos(x2)*exp(-a) +...
            -d*(x2/(pi*sqrt(x1^2+x2^2))));
    end
end

end

