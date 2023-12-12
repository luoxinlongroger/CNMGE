% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = hansen_function(x);
%
% This is a test problem subroutine used to test the performance of CNMGE.
% The test problem comes from [1].
%
% Input:
% x: It denotes the current x.
%
% Output: 
% func: It denotes the function value at x.
%
% grad_func: It denotes the gradient of the objective function at x. 
%
% References:
% [1] Andrei, N.: An unconstrained optimization test functions collection,
% Advanced Modeling and Optimization 10 (2008), 147-161.
%
function [func,grad_func] = hansen_function(x)
temp1 = 0;
temp2 = 0;
for i = 1 : 5
    temp1 = temp1 + i*cos((i+1)*x(2)+i);
    temp2 = temp2 + i*cos((i-1)*x(1)+i);
end
func = temp1*temp2;

% Compute the gradient of the objective function at x. 
if nargout > 1
    n = length(x);
    grad_func = zeros(n,1);
    temp1 = 0;
    temp2 = 0;
    temp3 = 0;
    temp4 = 0;
    for i = 1 : 5
        temp1 = temp1 - i*(i-1)*sin((i-1)*x(1)+i);
        temp2 = temp2 + i*cos((i+1)*x(2)+i);
        temp3 = temp3 + i*cos((i-1)*x(1)+i);
        temp4 = temp4 - i*(i+1)*sin((i+1)*x(2)+i);
    end
    grad_func(1) = temp1*temp2;
    grad_func(2) = temp3*temp4;
end 

end

