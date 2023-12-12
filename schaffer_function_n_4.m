%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SCHAFFER FUNCTION N. 4
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
% function [func,grad_func] = schaffer_function_n_4(x);
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
%
function [func,grad_func] = schaffer_function_n_4(x)
x1 = x(1);
x2 = x(2);

temp1 = (cos(sin(abs(x1^2-x2^2))))^2 - 0.5;
temp2 = (1 + 0.001*(x1^2+x2^2))^2;

func = 0.5 + temp1/temp2;

% Compute the gradient of the objective function at x.
if nargout > 1
    x1 = x(1);
    x2 = x(2);
    dif12 = x1^2-x2^2;
    dif21 = x2^2-x1^2;
    sin12 = sin(dif12);
    sin21 = sin(dif21);
    cos12 = cos(dif12);
    cos21 = cos(dif21);
    cs12 = cos(sin12);
    cs21 = cos(sin21);
    ss12 = sin(sin12);
    ss21 = sin(sin21);
    sum12 = 1+0.001*(x1^2+x2^2);
    n = 2;
    grad_func = zeros(n,1);
    if(dif12>=0)
        grad_func(1) = -4*x1*cos12*cs12*ss12/(sum12^2)-...
            0.004*x1*(cs12^2-0.5)/(sum12^3);
        grad_func(2) = 4*x2*cos12*cs12*ss12/(sum12^2)-...
            0.004*x2*(cs12^2-0.5)/(sum12^3);
    else
        grad_func(1) = 4*x1*cos21*cs21*ss21/(sum12^2)-...
            0.004*x1*(cs21^2-0.5)/(sum12^3);
        grad_func(2) = -4*x2*cos21*cs21*ss21/(sum12^2)-...
            0.004*x2*(cs21^2-0.5)/(sum12^3);
    end
end

end
