function [success,Plist,flist,alfa,radius,active,changes] = merge(x,f,...
         alfa_ini,radius_ini,Plist,flist,alfa,radius,...
         active,suf_decrease,poll,changes)
%
% Purpose:
%
%    Function merge compares a new evaluated point with the current list of
%    points, deciding if it should be added to it and updating the list.
%
% Input:  
%
%         x (Point to be compared.)
%
%         f (Corresponding objective function value.)
%
%         alfa_ini (Initial step size.)
%
%         radius_ini (Initial radius of comparison.)
%
%         Plist (Current list of points.)
%
%         flist (Corresponding objective function values.)
%
%         alfa (Corresponding step sizes.)
%
%         radius (Corresponding radius of comparison.)
%
%         active (Corresponding point status.) 
%
%         suf_decrease (0-1 variable: 1 if the algorithm uses a
%                      globalization strategy based in imposing a 
%                      sufficient decrease condition; 0 otherwise.)
% 
%         poll (0-1 variable: 1 if merging is performed inside the poll
%              step, 0 otherwise.)
%
%         changes (Record of points under analysis.)
%
% Output: 
%
%         success (0-1 variable: 1 if a better point was found; 0 otherwise.)
%
%         Plist (Updated list of points.)
%
%         flist (Corresponding objective function values.)
%
%         alfa (Corresponding step sizes.)
%
%         radius (Corresponding radius of comparison.)
%
%         active (Corresponding point status.)
%
%         changes (Record of points analysed.)
%
% Functions called: forcing (Provided by the optimizer.)
%
% GLODS Version 0.2
%
% Copyright (C) 2014 A. L. Custódio and J. F. A. Madeira.
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
success   = 0;
dist_list = sqrt(sum((Plist-repmat(x,1,size(Plist,2))).^2,1));
if min(dist_list-radius)>0
   success = 1;
   Plist   = [Plist,x];
   flist   = [flist,f];
   alfa    = [alfa,alfa_ini];
   radius  = [radius,radius_ini];
   active  = [active,1];
   changes = [changes,1];  
else
    if min(dist_list) ~= 0
       index      = find(dist_list-radius<=0);
       m_index    = size(index,2);
       active_new = 0;
       alfa_new   = 0;
       radius_new = 0;
       idom       = 0;
       pdom       = 0;
       icomp      = 0;
       for i=1:m_index
           if f < flist(index(i)) - forcing(alfa(index(i)),suf_decrease)
              icomp            = 1;
              idom             = idom + active(index(i));
              active(index(i)) = 0;
              if alfa(index(i)) > alfa_new
                 alfa_new   = alfa(index(i));
                 radius_new = radius(index(i));
              end   
           else
              if flist(index(i)) <= f - forcing(alfa(index(i)),suf_decrease)
                 pdom = 1;
              end
           end
       end
       if pdom == 0
           active_new = 1;
           success    = 1;
       end
       if (idom > 0) || (suf_decrease == 0 && pdom == 0 && icomp == 1)
          Plist   = [Plist,x];
          flist   = [flist,f];
          active  = [active,active_new];
          changes = [changes,1];
          if poll
             alfa   = [alfa,alfa(1)];
             radius = [radius,alfa(1)];  
          else
             alfa   = [alfa,alfa_new];
             radius = [radius,radius_new];    
          end
       end
    end
end
%
% End of merge.