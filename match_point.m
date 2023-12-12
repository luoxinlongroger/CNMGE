function [match,x,f] = match_point(x,xnorm,CacheP,CacheF,CachenormP,tol_match);
%
% Purpose:
%
%    Function match_point scans a list of previously evaluated points,
%    trying to match a point provided by the optimizer. When matching 
%    is successful, the function returns the corresponding objective
%    function value.
%
% Input: 
%
%         x (Point to be checked.)
%
%         xnorm (1-norm of the point to be checked.)
%
%         CacheP (Matrix of points to be used in comparisons,
%                storing the points columnwise.)
%
%         CacheF (Matrix storing the objective function values
%                corresponding to the points in CacheP.)
%
%         CachenormP (Vector storing the 1-norm of the points in CacheP.)
%
%         tol_match (Tolerance value within which two points are
%                   considered as equal.)
%
% Output: 
%
%         match (0-1 variable: 1 if x was previously evaluated; 0
%         otherwise.)
%
%         x (Vector storing the matched point coordinates.)
%
%         f (Objective function value at the matched point.)
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
% Prune the points that do not satisfy a 1-norm criterion.
%
f     = [];
index = find(abs(CachenormP - xnorm) <= tol_match);
if ~isempty(index)
   CacheP = CacheP(:,index);
   CacheF = CacheF(:,index);
%
% Finish search.
%
   nCacheP = size(CacheP,2);
   index   = find(max(abs(CacheP-repmat(x,1,nCacheP)),[],1) <= tol_match);
end
match = ~isempty(index);
%
% Retrieve the point coordinates and the objective function values. 
%
if match
   x = CacheP(:,index(1));
   f = CacheF(:,index(1));
end
%
% End of match_point.