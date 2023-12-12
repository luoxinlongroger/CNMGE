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
% #   See also Buckley#46 (n = 25) (p. 45).
% #   SIF input: Ph. Toint, Dec 1989.
%
% #   classification SUR2-AN-V-0
%
function func = chnrosnb(x)
%chnrosnb   12878663
n = 50;

alph = [1.25 1.40 2.40 1.40 1.75 1.20 ...
    2.25 1.20 1.00 1.10 1.50 1.60 ...
    1.25 1.25 1.20 1.20 1.40 0.50 ...
    0.50 1.25 1.80 0.75 1.25 1.40 ...
    1.60 2.00 1.00 1.60 1.25 2.75 ...
    1.25 1.25 1.25 3.00 1.50 2.00 ...
    1.25 1.40 1.80 1.50 2.20 1.40 ...
    1.50 1.25 2.00 1.50 1.25 1.40 ...
    0.60 1.50]';

term1 = x(1:n-1);
term2 = x(2:n);
term3 = alph(2:n);

func = sum (...
    ((term1-(term2.^2)).^2)*16.*(term3.^2) + ...
    ((term2-1.0).^2) );



end

