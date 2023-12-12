%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% HOLDER TABLE FUNCTION
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
% xx = [x1, x2]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = holder_table_function(x);
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
%
% func: It denotes the objective function. 
%
% grad_func: It denotes the gradient of the objective function at x. 
%
% References:
% [1] Surjanovic, S., Bingham, D.: Virtual library of simulation
% experiments: test functions and datasets, available at
% http://www.sfu.ca/~ssurjano, January 2020.
function [func,grad_func] = holder_table_function(x)
x1 = x(1);
x2 = x(2);

temp1 = sin(x1)*cos(x2);
temp2 = exp(abs(1 - sqrt(x1^2+x2^2)/pi));

func = -abs(temp1*temp2);

% Compute the gradient of the objective function at x. 
if nargout > 1
    x1 = x(1);
    x2 = x(2);
    sinx1 = sin(x1);
    sinx2 = sin(x2);
    cosx1 = cos(x1);
    cosx2 = cos(x2);
    sum12 = x1^2+x2^2;
    e1 = exp(1-sqrt(sum12)/pi);
    e2 = exp(sqrt(sum12)/pi-1);
    flag1 = sinx1*cosx2;
    flag2 = 1-sqrt(sum12)/pi;
    n = 2;
    grad_func = zeros(n,1);
    if(flag1>=0 && flag2>=0)
        grad_func(1) = -cosx1*cosx2*e1+sinx1*cosx2*(x1/(pi*sqrt(sum12)))*e1;
        grad_func(2) = sinx1*sinx2*e1+sinx1*cosx2*(x2/(pi*sqrt(sum12)))*e1;
    elseif(flag1>=0 && flag2<0)
        grad_func(1) = -cosx1*cosx2*e2-sinx1*cosx2*(x1/(pi*sqrt(sum12)))*e2;
        grad_func(2) = sinx1*sinx2*e2-sinx1*cosx2*(x2/(pi*sqrt(sum12)))*e2;
    elseif(flag1<0 && flag2>=0)
        grad_func(1) = cosx1*cosx2*e1-sinx1*cosx2*(x1/(pi*sqrt(sum12)))*e1;
        grad_func(2) = -sinx1*sinx2*e1-sinx1*cosx2*(x2/(pi*sqrt(sum12)))*e1;
    elseif(flag1<0 && flag2<0)
        grad_func(1) = cosx1*cosx2*e2+sinx1*cosx2*(x1/(pi*sqrt(sum12)))*e2;
        grad_func(2) = -sinx1*sinx2*e2+sinx1*cosx2*(x2/(pi*sqrt(sum12)))*e2;
    end
end

end

