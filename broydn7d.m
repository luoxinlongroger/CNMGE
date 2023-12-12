% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function [func,grad_func] = ackley_function(x);
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
% grad_func: It denotes the gradient of the objective function  at x.
%
% References:
% #   Source:
% #   Ph.L. Toint,
% #   "Some numerical results using a sparse matrix updating formula in
% #   unconstrained optimization",
% #   Mathematics of Computation, vol. 32(114), pp. 839-852, 1978.
%
% #   See also Buckley#84
% #   SIF input: Ph. Toint, Dec 1989.
%
% #   classification OUR2-AN-V-0
%
function func = broydn7d(x)
%broydn7d  12872726
n = 1000;

term1 = x(2:n-1);
term2 = x(1:n-2);
term3 = x(3:n);
term4 = x(1:(n/2));
term5 = x((n/2+1):n);

termVer1 = (-2*x(2)+1+(3-2*x(1))*x(1));
termVer2 = (1-term2 - 2*term3 + (3-2*term1).*term1);
termVer3 = (-x(n-1)+1+ (3-2*x(n)) *x(n));
termVer4 = (term4 + term5);


func = (sqrt(termVer1^2))^(7/3) ...
    + sum( (sqrt(termVer2.^2)).^(7/3) )...
    + (sqrt(termVer3^2))^(7/3)...
    + sum( (sqrt(termVer4.^2)).^(7/3) );


end

