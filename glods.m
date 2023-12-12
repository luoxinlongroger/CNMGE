function [glods_profile,Plist,flist,alfa,radius,func_eval,outmin] = glods(func_f,file_ini,x_ini,lbound,ubound)
%
% Purpose:
%
% Function glods applies a global derivative-free optimization method to
% solve the bound constrained problem:
%
%    min f(x) s.t. lbound <= x <= ubound,
%
% where x is a real vector of dimension n. The derivatives of the function
% f are not used, being provided only function values. 
%
% The user must provide: func_f (for f function values).
%
% Input:
%
%         func_f (Name of the file defining the objective function.)
%
%         file_ini (Name of the file used for initialization.)
%
%         x_ini (Point to be considered as initialization; only required 
%               when a file is not used for initialization, list is set 
%               equal to 0 and bounds are not provided.)
%
%         lbound (Lower bounds on the problem variables. Also used for
%                list initialization.)
%
%         ubound (Upper bounds on the problem variables. Also used for
%                list initialization.)
%
% Output:
%
%         glods_profile (Record of the evolution in the function value
%                       as function of the number of evaluations.)
%
%         Plist (List of current approximations to local minimizers.)
%
%         flist (List of function values corresponding to the current 
%               approximations to local minimizers.)
%
%         alfa (List of the corresponding step size parameters.)
%
%         radius (List of the corresponding comparison radius.)
%
%         func_eval (Total number of function evaluations.)
%
%         If in parameters_glods.m file output is not set equal to 0, then 
%         a text file, glods_report.txt, is stored at the current directory, 
%         which records the iteration results.
%
% Functions called: func_f (Application, user provided.)
%
%                   match_point, merge, parameters_glods, 
%                   search_step (Provided by the optimizer).
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
time = clock;
format long e;
warning off all;
%
% Load algorithmic strategies, parameters, constants, and tolerances
% values.
%
parameters_glods;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization Step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Read the file used for initialization.
%
if (list == 4)
   fpoints = fopen(file_ini,'r');
   aux     = str2num(fgetl(fpoints));
   n       = aux(1);
   m       = aux(2);
   aux     = str2num(fgetl(fpoints));
   for i = 1:n
      aux = str2num(fgetl(fpoints));
      for j = 1:m
         x_ini(i,j) = aux(j);
      end
   end
   aux = str2num(fgetl(fpoints));
   aux = str2num(fgetl(fpoints));
   for j = 1:m
      f_ini(j) = aux(j);
   end
   aux = str2num(fgetl(fpoints));
   aux = str2num(fgetl(fpoints));
   for j = 1:m
      alfa_list(j) = aux(j);
   end
   aux = str2num(fgetl(fpoints));
   aux = str2num(fgetl(fpoints));
   for j = 1:m
      radius_list(j) = aux(j);
   end
   fclose(fpoints);
end
%
% Define the problem size and the initial list of points.
%
if isempty(lbound) || isempty(ubound) || ((sum(isfinite(lbound),1)+...
   sum(isfinite(ubound),1))~= (2*size(lbound,1)))
   fprintf('Error: Initial complete variable bounds should be provided.\n\n');
   return
else
   n = size(lbound,1);
   if (list == 0) || (list == 4)
      if ~isempty(x_ini)
         Pini = [x_ini];
      else
         Pini = [(lbound + ubound)/2];
      end
   else
      if (user_list_size == 0)
         nPini = n;
      end
%   
      if list == 3
          if (nPini == 1)
             Pini = [(lbound + ubound)/2];
          else
             Pini = repmat(lbound,1,nPini)+repmat([0:(nPini-1)]/(nPini-1),n,1)...
                    .*repmat((ubound-lbound),1,nPini);
             Pini = [Pini,(lbound + ubound)/2];
          end
      end
   end 
end
%
func_eval = 0;
flist     = [];
match     = 0;
if cache ~= 0
   CacheP     = [];
   CachenormP = [];
   CacheF     = [];
end
%
% Set the seed for random strategies used in the initialization of 
% the list of points.
%
if (list == 1) || (list == 2)
   rand('state', sum(100*clock));
end
%
% Evaluate the initial iterate list.
%
%    
while isempty(flist)
   if (list == 1)
      Pini = [repmat(lbound,1,nPini)+lhsdesign(nPini,n)'.*...
             repmat((ubound-lbound),1,nPini)];
   else
      if (list == 2)
         Pini = [repmat(lbound,1,nPini)+rand(n,nPini).*...
             repmat((ubound-lbound),1,nPini)];
      end
   end
   for i = 1:size(Pini,2)   
      x_ini = Pini(:,i);
%
%  Check feasibility.
%   
      feasible = 1;
      bound    = [x_ini - ubound ; - x_ini + lbound];
      if sum(bound <= 0) ~= 2*n
         feasible = 0;
      end
%
      if feasible
          if (list~=4)
             if cache ~= 0
%                 
%     Check if the point was already evaluated.
%
                x_ininorm = norm(x_ini,1);
                if ~isempty(CacheP)
                   [match,x_ini,ftemp] = match_point(x_ini,x_ininorm,CacheP,...
                                         CacheF,CachenormP,tol_match);
                end
             end
             if ~match
%
%      Evaluate the point and store the corresponding values.
%
                ftemp     = feval(func_f,x_ini);
                func_eval = func_eval + 1;
                if cache ~= 0
                   CacheP     = [CacheP,x_ini];
                   CachenormP = [CachenormP,x_ininorm];
                   CacheF     = [CacheF,ftemp];
                end
             end
          else
             ftemp = f_ini(:,i);
          end
          if isfinite(ftemp)
             if isempty(flist)
                 flist = [ftemp];
                 Plist = [x_ini];
                 if (list ~= 4)
                     alfa   = [alfa_ini];
                     radius = [radius_ini];
                 else
                     alfa   = [alfa_list(i)];
                     radius = [radius_list(i)];
                 end
                 active = [1];
                 glods_profile(func_eval) = min(flist);
             else
                 if (list == 4)
                    alfa_aux   = alfa_list(i);
                    radius_aux = radius_list(i); 
                 else
                    alfa_aux   = alfa_ini;
                    radius_aux = radius_ini;
                 end
                 [~,Plist,flist,alfa,radius,active,~] = merge(x_ini,ftemp,...
                       alfa_aux,radius_aux,Plist,flist,alfa,radius,active,...
                       suf_decrease,0,[]);
                 glods_profile(func_eval) = min(flist);
             end
          end
      end
   end
%
% Check if the iterate list is not empty.
%
   if isempty(flist) && (list~=1) && (list~=2)
      fprintf('Error: The optimizer did not generate a feasible point\n');
      fprintf('or the initial point provided is not feasible.\n');
      fprintf('Please try list=1 or list=2 in parameters file.\n\n');
      return
   end
   if isempty(flist) && stop_feval && (func_eval >= max_fevals)
      fprintf('Error: The optimizer did not generate a feasible point,\n');
      fprintf('considering the budget of functions evaluations provided.\n\n');
      return
   end
end
%
% Set the seed for random generation of poll directions.
%
if (dir_dense == 1)
   rand('state',1);
end
%
if search_size == 0
    search_size = n;
end
halt            = 0;
iter            = 0;
iter_suc        = 0;
unsuc_consec    = 0;
grid_size       = 1;
label_grid_size = 1;
%
% Print the iteration report header.
%
if output     
   fprintf('Iteration Report: \n\n');
   fprintf('| iter  | success | #active points |    min fvalue   |    min alpha    |    max alpha    |\n');
   print_format = ['| %5d |   %2s    |    %5d       | %+13.8e | %+13.8e | %+13.8e |\n'];
   fprintf(print_format, iter, '--', sum(active), min(flist(logical(active))),min(alfa(logical(active))), max(alfa(logical(active))));
   fresult = fopen('glods_report.txt','w');   
   fprintf(fresult,'Iteration Report: \n\n');
   fprintf(fresult,'| iter  | success | #active points |    min fvalue   |    min alpha    |    max alpha    |\n');
   print_format = ['| %5d |   %2s    |    %5d       | %+13.8e | %+13.8e | %+13.8e |\n'];
   fprintf(fresult,print_format, iter, '--', sum(active), min(flist(logical(active))),min(alfa(logical(active))), max(alfa(logical(active))));
end
%
while (~halt)
%      
   func_iter   = 0;
   aux_success = 0;
   success     = 0;
   poll        = 1;
   search      = 0;
   changes     = zeros(1,size(flist,2)); 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Search Step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Check if the search step should be performed.
%
   if iter ~= 0 && search_option ~= 0 
      if search_freq_type == 0
        if search_freq == 0
           search = 1;
        else
           if unsuc_consec == search_freq
              search = 1;
           end
        end
      else
         index = find(alfa >= tol_active_points);
         aux_active = active(index);
         index = find(aux_active);
         if size(index,2) <= min_active_points
            search = 1;
         end
      end
   end
%
% Generate the new set of points to be evaluated.
%
   finite = 0;
   if search && ~halt
      while ~finite
         unsuc_consec = 0;
         [Psearch,grid_size,label_grid_size] = search_step(search_option,...
                                            search_size,lbound,ubound,...
                                            grid_size,label_grid_size);
         if ~isempty(Psearch)
             for i = 1:size(Psearch,2)
                xtemp = Psearch(:,i);         
%
%         Check feasibility.
%
                feasible = 1;
                bound    = [xtemp - ubound ; - xtemp + lbound];
                if sum(bound <= 0) ~= 2*n
                  feasible = 0;
                end
%
                if feasible
                   if cache ~= 0
%                 
%     Check if the point was already evaluated.
%
                      x_tempnorm = norm(xtemp,1);
                      [match,xtemp,ftemp] = match_point(xtemp,x_tempnorm,CacheP,...
                                            CacheF,CachenormP,tol_match);
                   end
                   if ~match  
%
%            Evaluate the point and store the corresponding values.
%
                      ftemp     = feval(func_f,xtemp);
                      func_eval = func_eval + 1;
                      func_iter = func_iter + 1;
                      if cache ~= 0
                         CacheP     = [CacheP,xtemp];
                         CachenormP = [CachenormP,x_tempnorm];
                         CacheF     = [CacheF,ftemp];
                      end
                   end
                   if isfinite(ftemp)
                      finite = 1;
                      [success,Plist,flist,alfa,radius,active,changes] = ...
                       merge(xtemp,ftemp,alfa_ini,radius_ini,Plist,flist,...
                       alfa,radius,active,suf_decrease,0,changes);
                      aux_success = aux_success + success;
                      glods_profile(func_eval) = min(flist);
                   end
                end
             end
             if stop_feval && (func_eval >= max_fevals)
                 halt = 1;
             end            
             if aux_success >0
                success = 1;
                poll    = 0;
             else
                success = 0;
                poll    = 1;
             end
         end
      end
   end
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Poll Step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
   if poll && ~halt
%        
%     Generate the positive basis.
%
      if (~dir_dense)
         D = [eye(n) -eye(n)];
      else
         v     = 2*rand(n,1)-1;
         [Q,R] = qr(v);
         if ( R(1) > 1 )
            D = Q * [ eye(n) -eye(n) ];
         else
            D = Q * [ -eye(n) eye(n) ];
         end
      end
%      
%  Reorder the list of points for the new iteration.
%         
     [flist,index] = sort(flist,'ascend');
     Plist         = Plist(:,index);
     alfa          = alfa(:,index);
     radius        = radius(:,index);
     active        = active(:,index);
%
    
     index1 = find(alfa >= tol_stop);
     index2 = find(alfa < tol_stop);
     Plist  = [Plist(:,index1),Plist(:,index2)];
     flist  = [flist(index1),flist(index2)];
     alfa   = [alfa(index1),alfa(index2)];
     radius = [radius(index1),radius(index2)];
     active = [active(index1),active(index2)];
%
     aux_active = active;
     Plist      = [Plist(:,logical(active)),Plist(:,~logical(active))];
     flist      = [flist(logical(active)),flist(~logical(active))];
     alfa       = [alfa(logical(active)),alfa(~logical(active))];
     radius     = [radius(logical(active)),radius(~logical(active))];
     active     = [aux_active(logical(active)),aux_active(~logical(active))]; 
% 
%
%     Poll using the positive basis.
%
      nd      = size(D,2);
      count_d = 1;
      changes = [1,zeros(1,size(flist,2)-1)];
      while ~success && (count_d <= nd)
         xtemp = Plist(:,1) + alfa(1) * D(:,count_d);
%
%        Check feasibility.
%
         feasible = 1;
         bound    = [xtemp - ubound ; - xtemp + lbound];
         if sum(bound <= 0) ~= 2*n
           feasible = 0;
         end
%
         if feasible  
            if cache ~= 0
%                 
%     Check if the point was already evaluated.
%
               x_tempnorm = norm(xtemp,1);
               [match,xtemp,ftemp] = match_point(xtemp,x_tempnorm,CacheP,...
                                     CacheF,CachenormP,tol_match);
            end
            if ~match 
%
%           Evaluate the point and store the corresponding values.
%
               ftemp     = feval(func_f,xtemp);
               func_eval = func_eval + 1;               
               if cache ~= 0
                  CacheP     = [CacheP,xtemp];
                  CachenormP = [CachenormP,x_tempnorm];
                  CacheF     = [CacheF,ftemp];
               end
            end
            if isfinite(ftemp)
               [success,Plist,flist,alfa,radius,active,changes] =...
                   merge(xtemp,ftemp,alfa_ini,radius_ini,Plist,flist,alfa,...
                   radius,active,suf_decrease,poll,changes);
               glods_profile(func_eval) = min(flist);
            end
         end
         count_d = count_d + 1;
      end
   end
% 
   if success
%       
%  Update the counter for successful iterations.
%          
      iter_suc     = iter_suc + 1;
      unsuc_consec = 0;
%
%  Update the step size parameter.
% 
      alfa(logical(changes+active == 2))   = alfa(logical(changes+active == 2))*gamma_par;
      radius(logical(changes+active == 2)) = max(radius(logical(changes+active == 2)),alfa(logical(changes+active == 2)));
   else
      unsuc_consec = unsuc_consec + 1; 
      alfa(logical(changes+active == 2)) = alfa(logical(changes+active == 2))*beta_par;
   end     
%
%  Check if the stopping criteria are satisfied.
%
   if stop_alfa && (sum(alfa(logical(active)) >= tol_stop)==0)
      halt = 1;
   end
   if stop_feval && (func_eval >= max_fevals)
      halt = 1;
   end
%
   iter = iter + 1;
%
%
%  Print the iteration report.
%
   if output 
     print_format = ['| %5d |    %1d    |    %5d       | %+13.8e | %+13.8e | %+13.8e |\n'];
     fprintf(print_format,iter,success,sum(active),min(flist(logical(active))),min(alfa(logical(active))),max(alfa(logical(active))));
     fprintf(fresult,print_format,iter,success,sum(active),min(flist(logical(active))),min(alfa(logical(active))),max(alfa(logical(active))));
   end
%
%  Print the current active points.
%
   if (output == 2)
      fglods = fopen('glods_partial_results.txt','w');
      n = size(Plist,1);
      m = sum(active);
      print_format = '%d %d\n\n';
      fprintf(fglods,print_format,n,m);
      print_format = [];   
      for i = 1:m
        print_format  = strcat(print_format,' %+21.16e ');
      end
      print_format  = strcat(print_format,'\n');
      for j=1:n
         fprintf(fglods,print_format,Plist(j,logical(active)));
      end
      fprintf(fglods,'\n');
      fprintf(fglods,print_format,flist(logical(active)));
      fprintf(fglods,'\n');
      fprintf(fglods,print_format,alfa(logical(active)));
      fprintf(fglods,'\n');
      fprintf(fglods,print_format,radius(logical(active)));
      fclose(fglods);
   end    
end
Plist = Plist(:,logical(active));
flist = flist(logical(active));
alfa  = alfa(logical(active));
radius = radius(logical(active));
time = etime(clock,time);
%
% Print the final report in the screen.
%
fprintf('\n Final Report: \n\n');
print_format = 'Elapsed Time = %10.3e \n\n';
fprintf(print_format,time);
fprintf('| #iter | #isuc | #active points | #fevals |    min fvalue   |\n');
print_format = ['| %5d | %5d |    %5d       |  %5d  | %+13.8e |\n\n'];
fprintf(print_format,iter,iter_suc,sum(active),func_eval,min(flist));
%
if output
%   
%  Print the final report in the results file.
%
   fprintf(fresult,'\n Final Report: \n\n');
   print_format = 'Elapsed Time = %10.3e \n\n';
   fprintf(fresult,print_format,time);
   fprintf(fresult,'| #iter | #isuc | #active points | #fevals |   min fvalue    |\n');
   print_format = ['| %5d | %5d |    %5d       |  %5d  | %+13.8e |\n\n'];
   fprintf(fresult,print_format,iter,iter_suc,sum(active),func_eval,min(flist));
   fclose(fresult);
   outmin = min(flist);
end


%
% End of glods.