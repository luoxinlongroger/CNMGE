% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = molecular_energy_problem(x);
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
% [1] Lavor, C., Maculan, N.: A function to test methods applied to global 
% minimization of potential energy of molecules, Numer. Algorithms 35
% (2004), 287-300.
%
function [func,grad_func] = molecular_energy_problem(x)
n = length(x);
k = 4.141720682;
r = 10.60099896;
p = 0;
q = 0;
for i = 1 : n
    p = p + 1 + cos(3*x(i));
    q = q + ((-1)^i)/sqrt(r-k*cos(x(i)));
end
func = p + q;

% Compute the gradient of the objective function at x. 
if nargout > 1
    k = 4.141720682;
    r = 10.60099896;
    n = length(x);
    grad_func = zeros(n,1);
    for i = 1 : n
      grad_func(i) = -3*sin(3*x(i))+...
        0.5*(-1)^(mod(i-1,2))*k*sin(x(i))/((r-k*cos(x(i)))*sqrt(r-k*cos(x(i))));
    end
end

end

