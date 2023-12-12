function [Psearch,grid_size,label_grid_size] = search_step(search_option,...
                    search_size,lbound,ubound,grid_size,label_grid_size)
%     
% Purpose:
%   
%    Function search_step generates a new set of points, to be used in the
%    search step, which should be asymptotically dense in the feasible region.
%    
% Input:  
%
%         search_option(0-5 variable: 0 if no search step should be performed;
%                      1 if a latin hypercube sampling strategy is considered; 
%                      2 if random sampling is used; 3 if points are uniformly 
%                      spaced in a grid (2^n-Centers strategy); 4 if points 
%                      are generated using Halton sequences; 5 if generation 
%                      is based on Sobol numbers.)
%
%         search_size (Number of points to be generated at each search
%                     step, when search_option is not set equal to 3. If 
%                     search_option is set equal to 3 then search_size is
%                     equal to 2^(problem dimension).) 
%
%         lbound (Lower bound on the problem variables.)
%
%         ubound (Upper bound on the problem variables.)
%
%         grid_size (Level of the current grid when search_option is set 
%                   equal to 3.)
%
%         label_grid_size (Label that indicates the current point to be
%                         used in the current grid when search_option is
%                         set equal to 3.)
%
% Output: 
%
%         Psearch (New list of points to be evaluated in the search step.)
%
%         grid_size (Level of the last used grid when search_option is set 
%                   equal to 3.)
%
%         label_grid_size (Label that indicates the next point to be used 
%                         in the grid when search_option is set equal to 3.)
%
% Functions called: replicate (Provided by the optimizer).
%
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
Psearch = [];
n       = size(lbound,1);
%
if (search_option == 1)
   Psearch = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
             lhsdesign(search_size,n)';
end
%
if (search_option == 2)
   Psearch = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
             rand(n,search_size);
end
%
if (search_option == 3)
    vaux       = [1/2^grid_size:1/2^(grid_size-1):1-1/2^grid_size]';
    Psearch    = vaux;
    for j = 1:n-1
       Psearch = replicate(vaux,Psearch);
    end
    if grid_size <= 2
       if grid_size == 1
          pindex = [1]; 
       else
          pindex = [1:2^n];
       end
       grid_size       = grid_size + 1;
       label_grid_size = 1;
    else
      index = [1]; 
      level = 0;
      count = size(index);
      while count < (2^n)^(grid_size-2)
         index_old = index;
         for i = 1:(2^(grid_size-2)-1)
           index_aux = (2^(grid_size-1))^level*i+index_old;
           index     = [index, index_aux];
         end
         count = size(index);
         level = level +1;
      end
      pindex = [index(label_grid_size) index(label_grid_size)+2^(grid_size-2)];
      level  = 1;
      while size(pindex) < 2^n  
        pindex = [pindex pindex+(2^(grid_size-1))^level*2^(grid_size-2)];
        level  = level +1;
      end
      if label_grid_size == (2^n)^(grid_size-2)
         grid_size       = grid_size + 1;
         label_grid_size = 1;
      else
         label_grid_size = label_grid_size + 1; 
      end
    end
    Psearch = Psearch(pindex,:)';
    if (grid_size == 2)
        vaux        = [0 1]';
        Psearch_aux = vaux;
        for j = 1:n-1
           Psearch_aux = replicate(vaux,Psearch_aux);
        end
        Psearch = [Psearch, Psearch_aux'];
    end
    Psearch = Psearch.*repmat((ubound-lbound),1,size(Psearch,2))+...
              repmat(lbound,1,size(Psearch,2));
end
%
if (search_option == 4)
    Lhalton   = haltonset(n,'Skip',n+1);
    Lhalton   = scramble(Lhalton,'RR2');
    Psearch   = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
                Lhalton(grid_size:grid_size + search_size - 1,:)';
    grid_size = grid_size + search_size;
end
%
if (search_option == 5)
    Lsobol    = sobolset(n,'Leap',2^n);
    Psearch   = repmat(lbound,1,search_size) + repmat((ubound-lbound),1,search_size).*...
                Lsobol(grid_size:grid_size + search_size - 1,:)';
    grid_size = grid_size + search_size;
end
%
% End of search_step.