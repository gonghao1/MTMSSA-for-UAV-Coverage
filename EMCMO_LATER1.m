function a=EMCMO_LATER1(dima,c2m,c3m,max_iter,iter_ratio,luratio,iter,N,Food_position,Salps_X,dim,sss,c1,c2,c3,wcc,dim_length,ub_ud,lb_ud,ub_total,lb_total,best_Solves1,fitness_set1,pos_xyz,dimv,dim_total,lx_set,ly_set,ux_set,uy_set,desxy,dimorder,Archive_size)
%  for i=1:N
% 
%         %index=0;
%         %neighbours_no=0;
%             dis_objective=[];
%             f_order=PMX_function(Salps_X{5,i},Food_position{5,1},c1);
%             %f_order=PMX_function(Food_position{5,1},Salps_X{5,i},c1);
%             %%
%             disarray_objective=[];
%             for g=1:length(f_order)-1,
%             dis_objective(g)=sign(pos_xyz(f_order(g+1),3)-pos_xyz(f_order(g),3));
%             disarray_objective=[disarray_objective ones(1,dimv).*dis_objective(g)];
%             end
% 
%         %% improved_section3
%         %if i<=N/2
%  if i==1
%             %disarray_obejective=[ones(1,1001).*dis_objective(i)];
%             %%  
%             %distance_result=calcu_dis(pos_xyz,f_order);
%             Salps_X{5,i}=f_order;
%             [lx_set,ux_set,ly_set,uy_set,desxy]=upper_lower_update(pos_xyz,Salps_X,dimv-1,N,luratio);% update upper and lower
%             for j=1:1:dim-1,% deal with 'dim' to 4
%                %c_total=[ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3]; 
%                c_total1=ones(1,dim_total(1)*10).*sss(j);
%                c_total2=ones(1,dim_total(1)*10).*c3;
%                temp_salps=[];
%                if j~=1 && j~=2
%                  for m=1:dim_length(j)*10,
%                 %if c3<0.5
%                 %% improved_section2
%                 if  disarray_objective(m)>0 && c_total2(m)>0.5
%                     temp_arr1=Food_position(j);
%                     temp_salps(m)=temp_arr1{1}(m)+c1*((ub_ud(1,j)-lb_ud(1,j))*c_total1(m)+lb_ud(1,j));
%                 elseif disarray_objective(m)>0 && c_total2(m)<0.5
%                     temp_arr1=Food_position(j);
%                     temp_salps(m)=temp_arr1{1}(m)-c1*((ub_ud(1,j)-lb_ud(1,j))*c_total1(m)+lb_ud(1,j));
%                 elseif disarray_objective(m)<0 && c_total2(m)>0.5
%                     temp_arr1=Food_position(j);
%                     temp_salps(m)=temp_arr1{1}(m)+c1*((ub_ud(2,j)-lb_ud(2,j))*c_total1(m)+lb_ud(2,j));
%                 elseif disarray_objective(m)<0 && c_total2(m)<0.5
%                     temp_arr1=Food_position(j);
%                     temp_salps(m)=temp_arr1{1}(m)-c1*((ub_ud(2,j)-lb_ud(2,j))*c_total1(m)+lb_ud(2,j));
%                 end
%                  end
%                elseif j==1 
%                    temp_arr1=Food_position(j);
%                   for m=1:dimorder,
%                      for n=1:dim_length(j),
%                      temp_salps(n+(dim_length(j)*(m-1)))=temp_arr1{1}(n+(dim_length(j)*(m-1)))+c1*((ux_set(i,m)-lx_set(i,m))*c_total1(n+(dim_length(j)*(m-1)))+lx_set(i,m));
%                      end
%                      for n=1:dim_length(j),
%                      temp_salps(n+(dim_length(j)*(m-1)))=temp_salps(n+(dim_length(j)*(m-1)))*desxy{i,1}(m)/sum(temp_salps(1+(dim_length(j)*(m-1)):dim_length(j)*m));
%                      end
%                   end 
%                elseif j==2
%                   temp_arr1=Food_position(j);
%                   for m=1:dimorder,
%                      for n=1:dim_length(j),
%                      temp_salps(n+(dim_length(j)*(m-1)))=temp_arr1{1}(n+(dim_length(j)*(m-1)))+c1*((uy_set(i,m)-ly_set(i,m))*c_total1(n+(dim_length(j)*(m-1)))+ly_set(i,m));
%                      end
%                      for n=1:dim_length(j),
%                      temp_salps(n+(dim_length(j)*(m-1)))=temp_salps(n+(dim_length(j)*(m-1)))*desxy{i,2}(m)/sum(temp_salps(1+(dim_length(j)*(m-1)):dim_length(j)*m));
%                      end
%                   end 
%                end
%                 %end
%                Salps_X{j,i}=temp_salps;
% 
%             end
%             %% adjustment order 
%             %% improved_section4
%             % Salps_X{5,i}=f_order;
%             % Salps_X{5,i}=PMX_function(Salps_X{5,i},Food_position{5,1},c1);
%         elseif i>=2 && i<N+1
%            % Salps_X{5,i}=PMX_function(Food_position{5,1},Salps_X{5,i},c1);
%             [lx_set,ux_set,ly_set,uy_set,desxy]=upper_lower_update(pos_xyz,Salps_X,dimv-1,N,luratio);% update upper and lower
%               temp_salps=[];
%                 if dominates(fitness_set1{iter}(i,:),fitness_set1{iter}(i-1,:))
%                     for v=1:dim-1,
%                         temp_salps=[];
% 
%               for m=1:dim_length(v)*dimorder,
% %             point1=Salps_X{v,i-1};
% %             point2=Salps_X{v,i};
% %             Salps_X{v,i}=(point2+point1)./2;
%             %% improved_section5
%               temp_arr1=best_Solves1{1,iter};
%               temp_salps(m)=(1-wcc)*temp_arr1{v,i}(m)+(wcc*temp_arr1{v,i-1}(m));
%               end
%                     if v~=1 && v~=2
%                         continue
%                     elseif v==1
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,1}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                     elseif v==2
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,2}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                      end
%               Salps_X{v,i}=temp_salps;
%                     end
%                 elseif dominates(fitness_set1{iter}(i-1,:),fitness_set1{iter}(i,:))
%                     for v=1:dim-1,
%                          temp_salps=[];
%               for m=1:dim_length(v)*dimorder,
%               temp_arr1=best_Solves1{1,iter};
%               temp_salps(m)=wcc*temp_arr1{v,i}(m)+((1-wcc)*temp_arr1{v,i-1}(m));
%               end
%                      if v~=1 && v~=2
%                         continue
%                     elseif v==1
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,1}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                     elseif v==2
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,2}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                      end
%                     Salps_X{v,i}=temp_salps;
%                     end
%                 else
%                     for v=1:dim-1,
%                          temp_salps=[];
%                 for m=1:dim_length(v)*dimorder,
%                 temp_arr1=best_Solves1{1,iter};
%                 temp_salps(m)=(1-wcc)*temp_arr1{v,i}(m)+(wcc*temp_arr1{v,i-1}(m));  
%                 end
%                      if v~=1 && v~=2
%                         continue
%                     elseif v==1
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,1}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                     elseif v==2
%                         for g=1:dimorder,
%                      for n=1:dim_length(v),
%                      temp_salps(n+(dim_length(v)*(g-1)))=temp_salps(n+(dim_length(v)*(g-1)))*desxy{i,2}(g)/sum(temp_salps(1+(dim_length(v)*(g-1)):dim_length(v)*g));
%                      end 
%                         end
%                      end
%                      Salps_X{v,i}=temp_salps;
%                     end
%                 end
% 
%             Salps_X{5,i}=f_order;
% end
% 
% %         for v=3:4,
% %         Flag4ub=Salps_X{v,i}>ub_total(v);
% %         Flag4lb=Salps_X{v,i}<lb_total(v);
% %         Salps_X{v,i}=(Salps_X{v,i}.*(~(Flag4ub+Flag4lb)))+ub_total(v).*Flag4ub+lb_total(v).*Flag4lb;  
% %         end
%         for v=3:4,
%             for l=1:length(dis_objective),
%                 if dis_objective(l)==1
%         Flag4ub=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)>ub_ud(1,v);
%         Flag4lb=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)<lb_ud(1,v);
%         Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)=(Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l).*(~(Flag4ub+Flag4lb)))+ub_ud(1,v).*Flag4ub+lb_ud(1,v).*Flag4lb;  
%                 else
%         Flag4ub=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)>ub_ud(2,v);
%         Flag4lb=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)<lb_ud(2,v);
%         Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)=(Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l).*(~(Flag4ub+Flag4lb)))+ub_ud(2,v).*Flag4ub+lb_ud(2,v).*Flag4lb;     
%                 end
%              end
%         end
% 
%            for v=1:2,
%             if v==1
%                 for l=1:dimorder,
%         Flag4ub=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))>ux_set(i,l);
%         Flag4lb=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))<lx_set(i,l);
%         Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*(~(Flag4ub+Flag4lb)))+ux_set(i,l).*Flag4ub+lx_set(i,l).*Flag4lb;  
%                 end
%             elseif v==2
%                 for l=1:dimorder,
%         Flag4ub=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))>uy_set(i,l);
%         Flag4lb=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))<ly_set(i,l);
%         Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*(~(Flag4ub+Flag4lb)))+uy_set(i,l).*Flag4ub+ly_set(i,l).*Flag4lb;  
%                 end
%             end
%             end
%         %Salps_tempfit(i,:)=ObjectiveFunction(Salps_X(:,i)',dim_length,penalty_constant*k,pos_xyz);
%   end
if iter<=max_iter*iter_ratio
[Salps_X, final_disobj]=mtmssa_update_motion2_update2(fitness_set1,dimv,sss,best_Solves1,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2,c3,dimorder,ub_ud,lb_ud,dim_length,dim_total,wcc,iter,Archive_size,max_iter);   
else
[Salps_X, final_disobj]=mssa_update_motion2(dim_total,dimv,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2m,c3m,dimorder,ub_ud,lb_ud,dim_length);
end
    a=Salps_X;
end