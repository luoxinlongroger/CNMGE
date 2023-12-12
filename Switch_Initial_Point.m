% This code is developed by Hang Xiao and Xinlong Luo. January, 2020.
% This code is revised by Hang Xiao and Xinlong Luo. December, 2021.
% luoxinlong@gmail.com or luoxinlong@bupt.edu.cn
% xiaohang0210@bupt.edu.cn
% The code comes with no guarantee or warranty of any kind.
%
% function x0_p = Switch_Initial_Point(x0, ini_number)
%
% This solver detetermines the global minmum of unconstrained optimization 
% problem by using continuation Newton methods with deflation techniques 
% and quasi-genetic ecolution. This program runs all programs successively 
% in desired way.
%
% This is the subroutine for switching initial points of this software package. 
%
% Have a look at title="Continuation Newton methods with deflation 
% techniques and quasi-genetic evolution for global optimization problems", 
% doi:10.21203/rs.3.rs-1102775/v1to, to see what this code does.
% 
% Input:
%
% x0: It denotes the initial point.
%
% ini_number: It denotes the serial number of the initial point.
%
% Output: 
% x0_p: It denotes the initial point selected for the continuation Newton method
% to find a stationary point of the objective function.
%
function x0_p = Switch_Initial_Point(x0, ini_number)
% xm is the length of vector x_k.
[xm,~] = size(x0);

% Initialize the next initial point.
x0_p = zeros(xm,1);

switch ini_number
    case 1
        % The first initial point is x0. 
        x0_p = x0;         
    case 2
        % The second initial point is -x0.
        x0_p = -x0;         
    case 3
        % The third initial point is [1,2,...,xm].
        ii = 1:xm;
        x0_p = ii';           
    case 4
        % The fourth initial point is [ones(n/2,1);-ones(n/2),1].
        if (xm == 1)
            x0_p(1) = 1;
        else 
            xmh = ceil(xm/2);
            x0_p(1:xmh) = 1;
            x0_p(xmh+1:xm) = -1;
        end         
    case 5
        % The fifth initial point is [-ones(n/2,1);ones(n/2),1].
        if (xm == 1)
            x0_p(1) = -1;
        else 
            xmh = ceil(xm/2);
            x0_p(1:xmh) = -1;
            x0_p(xmh+1:xm) = 1;
        end        
    case 6
        % The sixth initial point is [xm,xm-1,...,1];
        ii = xm:-1:1;
        x0_p = -ii'; 
%     case 7
%         % The seven initial point is [-1,-1,..., -1]. 
%         x0_p = -ones(xm,1);  
%     case 8
%         % The eigth initial point is [1,1,...,1]
%         x0_p = ones(xm,1);
end        
end

