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
% #   Source:  problem 32 in
% #   Ph. L. Toint,
% #   "Test problems for partially separable optimization and results
% #   for the routine PSPMIN",
% #   Report 83/4, Department of Mathematics, FUNDP (Namur, B), 1983.
%
% #   See  also Buckley#18
% #   SIF input: Ph. Toint, Dec 1989.
%
% #   classification OUR2-AY-V-0
%
function func = cragglvy(x)
%cragglvy  12878765(1000) 12878694(5000)
m = 499;
n = 2*m+2;

term1 = x(1:2:n-2); %2i-1
term2 = x(2:2:n-2); %2i
term3 = x(3:2:n);   %2i+1
term4 = x(4:2:n);   %2i+2

func = sum( ...
    ((exp(term1)-term2).^4) + ...
    100*((term2-term3).^6) + ...
    ((tan(term3-term4)+term3-term4).^4) + ...
    ((term1).^8) + ...
    ((term4-1.0).^2) );



end

