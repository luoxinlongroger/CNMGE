function [rho] = forcing(alfa,suf_decrease)
%
% Purpose:
%
%    Function forcing implements a globalization strategy based on
%    imposing a sufficient decrease condition.
%
% Input:  
%
%         alfa (Step size, given by the optimizer.)
%
%         suf_decrease (0-1 variable: 1 if the algorithm uses a
%                      globalization strategy based in imposing a 
%                      sufficient decrease condition; 0 otherwise.)
%
% Output: 
%
%         rho (Forcing function value at the given step size.)        
%
% GLODS Version 0.1
%
% Copyright (C) 2013 A. L. Custódio and J. F. A. Madeira.
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 3.0 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc.,51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
%
%
if suf_decrease
   rho = alfa^2;
else
   rho = 0;
end
%
% End of forcing.