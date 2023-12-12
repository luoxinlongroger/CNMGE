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
% #   L.C.W. Dixon and Z. Maany,
% #   "A family of test problems with sparse Hessians for unconstrained
% #   optimization",
% #   TR 206, Numerical Optimization Centre, Hatfield Polytechnic, 1988.
%
% #   See also Buckley#221 (p. 49)
% #   SIF input: Ph. Toint, Dec 1989.
% #              correction by Ph. Shott, January 1995.
%
% #   classification OUR2-AN-V-0
%
function func = dixmaank(x)
%dixmaank   其余的仅参数不同 12952012
m = 1000;
n = 3*m;

alpha = 1.0;
beta = 0.125;
gamma = 0.125;
delta  = 0.125;
K = [2 0 0 2]';

var1 = (1:n)';
var2 = (1:n-1)';
var3 = (1:2*m)';
var4 = (1:m)';

xTerm1 = x(1:n-1);  %1,n-1
xTerm2 = x(2:n);    %2,n
xTerm3 = x(1:2*m);
xTerm4 = x((1+m):(2*m+m));
xTerm5 = x(1:m);
xTerm6 = x((1+2*m):(m+2*m));

func =  1 + ...
    sum( alpha* ((x.^2).*((var1/n).^K(1))) ) + ...
    sum( beta* ((xTerm1.^2).*((xTerm2+(xTerm2.^2)).^2) .* ((var2/n).^K(2))) )+ ...
    sum( gamma* ((xTerm3.^2) .* (xTerm4.^4) .* ((var3/n).^K(3))) ) + ...
    sum( delta* (xTerm5.*xTerm6.*((var4/n).^K(4))) );

end

