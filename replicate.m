function  [M_new] = replicate(v,M)
%
% Purpose:
%
%    Auxiliary function that computes a new matrix by duplicating a 
%    provided matrix and duplicating a provided vector.
%
% Input:  
%
%         v (Vector to be duplicated.)
%
%         M (Matrix to be duplicated.)
%
% Output: 
%
%         M_new (New matrix.)        
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
M_new = [];
   for i = 1:size(v)
    M_new = [M_new; M repmat(v(i),size(M,1),1)];
   end
end
%
% End of replicate.