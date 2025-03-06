function Positions=initialization2(SearchAgents_no,dim,ub,lb,M,diml,ub_ud,lb_ud,pos_xyz)
%% 
% clear all
% close all;
% clc
% dim=5;dimv=501;dima=dimv-1;dimorder=10;% ？？？变量维度
% N=200;
% lb=0;% 变量下界
% ub=1;% 变量上界
% obj_no=2;%目标函数数量
% if size(ub,2)==1
%     ub=ones(1,dim)*ub;
%     lb=ones(1,dim)*lb;
% end
% lb_ah=-3;ub_ah=3;lb_av=-1;ub_av=1;lb_vh=0.01;ub_vh=20;lb_vv=-3;ub_vv=5;
% if size(lb_vh,2)==1
%     lb_ah=ones(1,dima).*lb_ah;
%     ub_ah=ones(1,dima).*ub_ah;
%     lb_av=ones(1,dima).*lb_av;
%     ub_av=ones(1,dima).*ub_av;
%     lb_vh=ones(1,dimv).*lb_vh;
%     ub_vh=ones(1,dimv).*ub_vh;
%     lb_vv=ones(1,dimv).*lb_vv;
%     ub_vv=ones(1,dimv).*ub_vv;
% end
% dim_total=[dimv dima];
% dim_length=[dima dima dimv dimv];
% dim_length2=[dima dima dimv dimv dimorder];
% ub_total=[ub_ah(1),ub_av(1),ub_vh(1),ub_vv(1)];
% lb_total=[lb_ah(1),lb_av(1),lb_vh(1),lb_vv(1)];
% ub_up=[3,1,20,5];
% lb_up=[-3,-1,0.01,0.01];
% ub_down=[3,1,20,0.01];
% lb_down=[-3,-1,0.01,-3];
% ub_ud=[ub_up;ub_down];
% lb_ud=[lb_up;lb_down];
% %% position xyz
% pos_x=2000*rand(1,dimorder+1);
% pos_y=2000*rand(1,dimorder+1);
% pos_z=10+100*rand(1,dimorder+1);
% pos_xyz=[pos_x;pos_y;pos_z]';
%%
Boundary_no=size(ub,2); % numnber of boundaries

% if Boundary_no==1
%     ub_new=ones(1,dim)*ub;
%     lb_new=ones(1,dim)*lb;
% else
%      ub_new=ub;
%      lb_new=lb;   
% end

% If each variable has a different lb and ub
%N,dim-1,ub_total,lb_total,dimorder+1,dim_length,ub_ud,lb_ud,pos_xyz


% SearchAgents_no=N;
% dim=dim-1;ub=ub_total;lb=lb_total;M=dimorder+1;diml=dim_length;

%xx1=order_select{1};
%xxx=1:1:M;
for i=1:SearchAgents_no,
    for j=1:diml(1)*30,
        order_select{j}=randperm(M);
        if sum(order_select{j})==sum(1:1:M)
        order_pos(j)=0;
        else 
            order_pos(j)=1;
        end
    end
%Positions{dim+1,i}=randperm(M);
findorder=find(order_pos==0);
Positions{dim+1,i}=order_select{findorder(1)};
[dista,distb]=calcu_dis(pos_xyz,Positions{dim+1,i});
distance{i}=[dista',distb'];
end
Orderrr=[];
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
       if i==dim-1
       for j=1:SearchAgents_no,
            temp_positions=[];dis_objective=[];
            for g=1:length(Positions{dim+1,j})-1,
            dis_objective(g)=sign(pos_xyz(Positions{dim+1,j}(g+1),3)-pos_xyz(Positions{dim+1,j}(g),3));
            end
               for g=1:length(dis_objective),
                  if dis_objective(g)>0
%                       [ts1,ts2]=sort(temp_salps(1+1001*(g-1):1001*g));
%                       temp_salps(1+1001*(g-1):1001*g)=ts1;
%                       ts3=Salps_X{2,i}(1+1001*(g-1):1001*g);% 
%                       Salps_X{2,i}(1+1001*(g-1):1001*g)=ts3(ts2);
                        randmm=weierstrass_function(rand(1,diml(i)),50);
                        randmm=abs(randmm)/2;
                        temp1=randmm.*(ub_ud(1,i)-lb_ud(1,i))+lb_ud(1,i);
                        %temp_positions(1+diml(i)*(g-1):diml(i)*g)=sort(temp1);
                        temp_positions(1+diml(i)*(g-1):diml(i)*g)=temp1;
                  else 
%                       [ts1,ts2]=sort(temp_salps(1+1001*(g-1):1001*g),'descend');
%                       temp_salps(1+1001*(g-1):1001*g)=ts1;
%                       ts3=Salps_X{2,i}(1+1001*(g-1):1001*g);
%                       Salps_X{2,i}(1+1001*(g-1):1001*g)=ts2(ts2);
                        randmm=weierstrass_function(rand(1,diml(i)),50);
                        randmm=abs(randmm)/2;
                        temp1=randmm.*(ub_ud(2,i)-lb_ud(2,i))+lb_ud(2,i);
                        %temp_positions(1+diml(i)*(g-1):diml(i)*g)=sort(temp1,'descend');  
                        temp_positions(1+diml(i)*(g-1):diml(i)*g)=temp1;  
                  end
               end
               Positions{i,j}=temp_positions;
               %temp_positions=[];dis_objective=[];
               Orderrr{1,j}=dis_objective;
       end
       else 
         for j=1:SearchAgents_no,
        temp_positions=rand(diml(i)*(M-1),SearchAgents_no).*(ub_i-lb_i)+lb_i;
        %Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
        Positions{i,j}=temp_positions(:,j)';
         end
        
       end

    end
   for i=1:size(Orderrr,2),
         t_array=[];
       for j=1:length(Orderrr{1,i})
       t_array=[t_array ones(1,diml(3)).*Orderrr{1,i}(j)];
       Orderrr{1,i}=t_array;
       end
   end
   %% 1
    % for i=1:size(Positions,2),
    % for j=1:length(Positions{1,1})/diml(1),
    % for k=1+(j-1)*diml(4):(diml(4)*j)-1,
    %     if Positions{4,i}(k+1)>Positions{4,i}(k)
    %        Positions{5,i}(k-(j-1))=abs(Positions{5,i}(k-(j-1)));
    %     elseif Positions{4,i}(k+1)<Positions{4,i}(k)
    %         Positions{5,i}(k-(j-1))=-abs(Positions{5,i}(k-(j-1)));
    %     else
    %         Positions{5,i}(k)=0;
    %     end
    %     % if Positions{4,i}(k+1)>Positions{4,i}(k) && Orderrr{1,i}(k)>0
    %     %    Positions{2,i}(k-(j-1))=abs(Positions{2,i}(k-(j-1)));
    %     % elseif Positions{4,i}(k+1)<Positions{4,i}(k) && Orderrr{1,i}(k)<0 
    %     %     Positions{1,i}(k-(j-1))=abs(Positions{1,i}(k-(j-1)));
    %     % elseif Positions{4,i}(k+1)>Positions{4,i}(k) && Orderrr{1,i}(k)<0
    %     %     Positions{2,i}(k-(j-1))=-abs(Positions{2,i}(k-(j-1)));
    %     % elseif Positions{4,i}(k+1)<Positions{4,i}(k) && Orderrr{1,i}(k)>0
    %     %     Positions{2,i}(k-(j-1))=-abs(Positions{2,i}(k-(j-1)));
    %     % elseif Positions{4,i}(k+1)==Positions{4,i}(k)
    %     %     Positions{2,i}(k-(j-1))=0;
    %     % end
    % end
    % end
    % end
    %% 2

%Positions=Positions';

end
