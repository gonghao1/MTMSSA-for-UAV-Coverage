function a=EMCMO_EARLY1(dima,c2m,c3m,max_iter,iter_ratio,luratio,iter,N,Food_position,Salps_X,dim,sss,c1,c2,c3,wcc,dim_length,ub_ud,lb_ud,ub_total,lb_total,best_Solves1,fitness_set1,pos_xyz,dimv,dim_total,lx_set,ly_set,ux_set,uy_set,desxy,dimorder,Archive_size)
%     for i=1:N
%             dis_objective=[];
%             %f_order=PMX_function(Salps_X{5,i},Food_position{5,1},c1);
%             f_order=PMX_function(Food_position{5,1},Salps_X{5,i},c1);
%             %%
%             disarray_objective=[];
%             for g=1:length(f_order)-1,
%             dis_objective(g)=sign(pos_xyz(f_order(g+1),3)-pos_xyz(f_order(g),3));
%             disarray_objective=[disarray_objective ones(1,dimv).*dis_objective(g)];
%             end
%  if i==1
% 
%             for j=1:1:dim-1,
%               %c_total=[ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3;ones(1,dim_total(1)*10).*c2;ones(1,dim_total(1)*10).*c3]; 
%                c_total1=ones(1,dim_total(1)*10).*sss(j);
%                c_total2=ones(1,dim_total(1)*10).*c3;
%                 temp_salps=[];
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
%                 end 
% 
%                 %end
%                Salps_X{j,i}=temp_salps;
%             end
%              Salps_X{5,i}=PMX_function(Salps_X{5,i},Food_position{5,1},c1);
%         elseif i>=2 && i<N+1
%               temp_salps=[];
%                 if dominates(fitness_set1{iter}(i,:),fitness_set1{iter}(i-1,:))
%                     for v=1:dim-1,
%                          temp_salps=[];
%               for m=1:dim_length(v)*10,
%               temp_arr1=best_Solves1{1,iter};
%               temp_salps(m)=(1-wcc)*temp_arr1{v,i}(m)+(wcc*temp_arr1{v,i-1}(m));
%               end
%               Salps_X{v,i}=temp_salps;
%                     end
%                 elseif dominates(fitness_set1{iter}(i-1,:),fitness_set1{iter}(i,:))
%                     for v=1:dim-1,
%                          temp_salps=[];
%               for m=1:dim_length(v)*10,
%               temp_arr1=best_Solves1{1,iter};
%               temp_salps(m)=wcc*temp_arr1{v,i}(m)+((1-wcc)*temp_arr1{v,i-1}(m));
%               end
%                     Salps_X{v,i}=temp_salps;
%                     end
%                 else
%                     for v=1:dim-1,
%                          temp_salps=[];
%                 for m=1:dim_length(v)*10,
%                 temp_arr1=best_Solves1{1,iter};
%                 temp_salps(m)=(1-wcc)*temp_arr1{v,i}(m)+(wcc*temp_arr1{v,i-1}(m));  
%                 end
%                      Salps_X{v,i}=temp_salps;
%                     end
%                 end
%                  
%             % Eq. (3.4) in the paper
%             %Salps_X{5,i}=randperm(dimorder+1);
%             Salps_X{5,i}=PMX_function(Salps_X{5,i},Food_position{5,1},c1);
% end
% % 
% %         Flag4ub=Salps_X{1:4,i}>ub_total';
% %         Flag4lb=Salps_X{1:4,i}<lb_total';
% %         Salps_X{1:4,i}=(Salps_X{1:4,i}.*(~(Flag4ub+Flag4lb)))+ub_total'.*Flag4ub+lb_total'.*Flag4lb;
% %         YUUU=[1,4,6,7]';YUUU_1=[4,3,2,1]';
% %         flag11=YUUU>1;
%         for v=1:4,
%         Flag4ub=Salps_X{v,i}>ub_total(v);
%         Flag4lb=Salps_X{v,i}<lb_total(v);
%         Salps_X{v,i}=(Salps_X{v,i}.*(~(Flag4ub+Flag4lb)))+ub_total(v).*Flag4ub+lb_total(v).*Flag4lb;  
%         end
%         %Salps_tempfit(i,:)=ObjectiveFunction(Salps_X(:,i)',dim_length,penalty_constant*k,pos_xyz);
%     end

% Salps_X2=EMCMO_EARLY(iter,N,Food_position2,Salps_X2,dim,sss,c1,c3,wcc,dim_length,ub_ud1,lb_ud1,ub_total1,lb_total1,best_Solves2,fitness_set2,pos_xyz,dimv,dim_total,lx_set1,ly_set1,ux_set1,uy_set1,desxy1,dimorder);
%Food_position=Food_position2;Salps_X=Salps_X2;ub_ud=ub_ud1;lb_ud-lb_ud1;ub_total=ub_total1;lb_total=lb_total1;best_Solves1=best_Solves2;fitness_set1=fitness_set2;lx_set=lx_set1;ly_set=ly_set1;ux_set=ux_set1;uy_set=uy_set1;desxy=desxy1;
           
%     fo
if iter<=max_iter*iter_ratio
[Salps_X, final_disobj]=mtmssa_update_motion2_update2(fitness_set1,dimv,sss,best_Solves1,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2,c3,dimorder,ub_ud,lb_ud,dim_length,dim_total,wcc,iter,Archive_size,max_iter);   
else
[Salps_X, final_disobj]=mssa_update_motion2(dim_total,dimv,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2m,c3m,dimorder,ub_ud,lb_ud,dim_length);
end
    a=Salps_X;
end