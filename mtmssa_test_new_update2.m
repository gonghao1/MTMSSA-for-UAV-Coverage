% clear
% close all
% clc
% Change these details with respect to your problem
%load('C:\Users\Administrator\Desktop\MSSA_MSSA_update\mssa1.mat');
function [fitness_set,best_Solves]=mtmssa_test_new_update2(betta11,iter_ratio,dimv,pos_xyz,max_iter),
ObjectiveFunction=@ZDT10_new; % 目标函数
ObjectiveFunction2=@ZDT10_new; % unconstraint function 
dim=5;
%dimv=5;
dima=dimv-1;dimorder=size(pos_xyz,1)-1;% ？？？变量维度
lb=0;% 变量下界
ub=1;% 变量上界
obj_no=2;%目标函数数量
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end
lb_ah=-1;ub_ah=1;lb_av=-1;ub_av=1;lb_vh=0.5;ub_vh=20;lb_vv=-3;ub_vv=5;
if size(lb_vh,2)==1
    lb_ah=ones(1,dima).*lb_ah;
    ub_ah=ones(1,dima).*ub_ah;
    lb_av=ones(1,dima).*lb_av;
    ub_av=ones(1,dima).*ub_av;
    lb_vh=ones(1,dimv).*lb_vh;
    ub_vh=ones(1,dimv).*ub_vh;
    lb_vv=ones(1,dimv).*lb_vv;
    ub_vv=ones(1,dimv).*ub_vv;
end
dim_total=[dimv dima];
dim_length=[dima dima dimv dimv];
dim_length2=[dimv dimv dimv dimv dimorder+1];
% ub_total=[ub_ah(1),ub_av(1),ub_vh(1),ub_vv(1)];
% lb_total=[lb_ah(1),lb_av(1),lb_vh(1),lb_vv(1)];
acceler_bound=[lb_ah(1),ub_ah(1),lb_av(1),ub_av(1)];
%% multitasking 
%ub_total1=[5,2,25,6];
% ub_total1=[1,1,40,7];
% lb_total1=[-5,-2,0.01,-5];
%lb_total1=[-5,-2,0.01,-3];
% ub_up=[3,1,0,sqrt(20^2+5^2),sqrt(1+1)];
% ub_up1=[5,2,0,sqrt(30^2+7^2),sqrt(1.5^2+1.5^2)];
% lb_up=[-3,-1,0,sqrt(0.5^2+0.1^2),-sqrt(1+1)];
% lb_up1=[-5,-2,0,sqrt(0.1^2+0.1^2),-sqrt(1.5^2+1.5^2)];
% ub_down=[3,1,0,sqrt(20^2+3^2),sqrt(1+1)];
% ub_down1=[5,2,0,sqrt(30^2+4^2),sqrt(1.5^2+1.5^2)];
% lb_down=[-3,-1,0,sqrt(0.5^2+0.05^2),-sqrt(1+1)];
% lb_down1=[-5,-2,0,sqrt(0.1^2+0.05^2),-sqrt(1.5^2+1.5^2)];

ub_up=[3,1,0,sqrt(20^2+5^2)];
ub_up1=[5,2,0,sqrt(35^2+7^2)];
lb_up=[-3,-1,0,sqrt(0.5^2+0.1^2)];
lb_up1=[-5,-2,0,sqrt(0.1^2+0.1^2)];
ub_down=[3,1,0,sqrt(20^2+3^2)];
ub_down1=[5,2,0,sqrt(40^2+4^2)];
lb_down=[-3,-1,0,sqrt(0.5^2+0.05^2)];
lb_down1=[-5,-2,0,sqrt(0.1^2+0.05^2)];

acceler_bound=[-sqrt(1+1),sqrt(1+1),-sqrt(1+1),sqrt(1+1)];
ub_ud=[ub_up;ub_down];
ub_ud1=[ub_up1;ub_down1];
lb_ud=[lb_up;lb_down];
lb_ud1=[lb_up1;lb_down1];

posx_bound=[];posy_bound=[];posz_bound=[];% 3-D position bound 
for t=1:size(pos_xyz,1)-1,
posx_bound(t,:)=[pos_xyz(t,1),pos_xyz(t+1,1)];
posy_bound(t,:)=[pos_xyz(t,2),pos_xyz(t+1,2)];
posz_bound(t,:)=[pos_xyz(t,3),pos_xyz(t+1,3)];
end

ub_total=[posx_bound(1,2),posy_bound(1,2),posx_bound(1,2),sqrt(20^2+5^2)];
lb_total=[posx_bound(1,2),posy_bound(1,2),posx_bound(1,2),sqrt(0.5^2+0.05^2)];
ub_total1=[posx_bound(1,2),posy_bound(1,2),posx_bound(1,2),sqrt(35^2+7^2)];
lb_total1=[posx_bound(1,2),posy_bound(1,2),posx_bound(1,2),sqrt(0.1^2+0.01^2)];
dim_length=[dimv dimv dimv dimv];
%dim_length2=[dimv dimv dimv dimv dimv dimorder];

%max_iter=300; % 迭代数量
N=100; % 种群数量
ArchiveMaxSize=500;
epslison=1;
%betta11=0.7;% multitasking
FES=0;MAXFES=max_iter*N;% multitasking
FES=FES+2*N;% multitasking
max_runtimes=1;
penalty_constant=2;
%iter_ratio=0.6;

c2_libm=zeros(1,max_iter);c2_libme=zeros(1,max_iter);c2_libml=zeros(1,max_iter);
c3_libm=zeros(1,max_iter);c3_libme=zeros(1,max_iter);c3_libml=zeros(1,max_iter);
c2_libm(1)=rand();c3_libm(1)=0.234563;c2_libme(1)=rand();c2_libml(1)=rand();c3_libme(1)=rand();c3_libml(1)=rand();
for i=1:max_iter-1,
c2_libm(i+1)=sin(pi*c2_libm(i));
c2_libme(i+1)=sin(pi*c2_libme(i));
c2_libml(i+1)=sin(pi*c2_libml(i));
if c3_libm(i)~=0
c3_libm(i+1)=mod((1/c3_libm(i)),1);
c3_libme(i+1)=mod((1/c3_libme(i)),1);
c3_libml(i+1)=mod((1/c3_libml(i)),1);
else 
    c3_libm(i+1)=1;
    c3_libme(i+1)=1;
    c3_libml(i+1)=1;
end
end
%% 
%% 
%% algorithm_begin
for k=1:max_runtimes,
%% archive
%Archive_X=zeros(100,dim);
for i=1:dim,
for j=1:ArchiveMaxSize,
Archive_X{i,j}=zeros(1,dim_length2(i));
Archive_X1{i,j}=zeros(1,dim_length2(i));
end
end
Archive_X=Archive_X';
Archive_X1=Archive_X1';
%Archive_X=zeros(100,dim);
%Archive_F=ones(100,obj_no)*inf;
Archive_F=ones(ArchiveMaxSize,obj_no)*inf;
Archive_F1=ones(ArchiveMaxSize,obj_no)*inf;
Archive_member_no=0;
Archive_member_no1=0;
%r=(ub-lb)/2;
%V_max=(ub(1)-lb(1))/10;
Food_fitness=inf*ones(1,obj_no);
Food_fitness2=inf*ones(1,obj_no);
Food_position=zeros(dim,1);
Food_position2=zeros(dim,1);
%Salps_X=initialization(N,dim,ub,lb);
Salps_X=initialization2(N,dim-1,ub_total,lb_total,dimorder+1,dim_length,ub_ud,lb_ud,pos_xyz);
Salps_X2=initialization2(N,dim-1,ub_total1,lb_total1,dimorder+1,dim_length,ub_ud1,lb_ud1,pos_xyz);% multitasking 
%fitness=zeros(N,2);
%V=initialization(N,dim,ub,lb);
%position_history=zeros(N,max_iter,dim);
%% 
%% improved_section6 更新c2和c3
c2_lib=zeros(1,max_iter);c2_libd=zeros(1,max_iter);c2_libe=zeros(1,max_iter);
c3_lib=zeros(1,max_iter);c3_libe=zeros(1,max_iter);c3_libl=zeros(1,max_iter);
c2_lib(1)=rand();c3_lib(1)=rand();c2_libd(1)=rand();c2_libe(1)=rand();c3_libe(1)=rand();c3_libl(1)=rand();
miu=70;
for i=1:max_iter-1,
%c2_lib(i+1)=sin(pi*c2_lib(i)); % sine 映射
c2_libd(i+1)=sin(miu*pi*c2_libd(i));
c2_libe(i+1)=sin(miu*pi*c2_libe(i));
c2_lib(i+1)=mod(c2_libd(i+1)+c2_libe(i+1),1);     % 改进sine 映射
c3_lib(i+1)=2.595*c3_lib(i)*(1-c3_lib(i)^2); % cubic 映射-> main
c3_libe(i+1)=2.595*c3_libe(i)*(1-c3_libe(i)^2); %  cubic映射-> early phase
c3_libl(i+1)=2.595*c3_libl(i)*(1-c3_libl(i)^2); %  cubic映射-> later phase
end
%%　improved_section7 levy flight mechanism
betaa=0.3+1.69*rand();% [0.3,1.99]
betaa1=0.3+1.69*rand();% early phase
betaa2=0.3+1.69*rand();% later phase
betaa=1.5;betaa1=1.5;betaa2=1.5;
%% +++++++++++++++++++++++++ initialization distance +++++++++++++++++++++++++++++++++
desxy={
    };
desxy1={
    };
n_gauss=12;% unused
for f=1:N,

  posxy_bound_temp=pos_xyz(Salps_X{dim,f},:);
  posxy_bound_temp1=pos_xyz(Salps_X2{dim,f},:);
  luratio=0.9; % adjustment parameter
  for d=1:size(pos_xyz,1)-1,
  lx_set(f,d)=posxy_bound_temp(d,1);% lower bound 
  ux_set(f,d)=posxy_bound_temp(d+1,1);% upper bound
  ly_set(f,d)=posxy_bound_temp(d,2);
  uy_set(f,d)=posxy_bound_temp(d+1,2);
  lz_set(f,d)=posxy_bound_temp(d,3);
  uz_set(f,d)=posxy_bound_temp(d+1,3);
  lx_set1(f,d)=posxy_bound_temp1(d,1);
  ux_set1(f,d)=posxy_bound_temp1(d+1,1);
  ly_set1(f,d)=posxy_bound_temp1(d,2);
  uy_set1(f,d)=posxy_bound_temp1(d+1,2);
  lz_set1(f,d)=posxy_bound_temp1(d,3);
  uz_set1(f,d)=posxy_bound_temp1(d+1,3);
  for r=1:dimv,
  %desp{f,1}(d)=(sqrt(12/n_gauss)*(sum(rand(1,n_gauss))-(n_gauss/2)))+(desx/(dimv-1));
  desp{f,1}(d,r)=lx_set(f,d)+(abs(weierstrass_function(rand(),50))*(ux_set(f,d)-lx_set(f,d)));
  desp{f,2}(d,r)=ly_set(f,d)+(abs(weierstrass_function(rand(),50))*(uy_set(f,d)-ly_set(f,d)));
  desp{f,3}(d,r)=lz_set(f,d)+(abs(weierstrass_function(rand(),50))*(uz_set(f,d)-lz_set(f,d)));
  desp1{f,1}(d,r)=lx_set1(f,d)+(abs(weierstrass_function(rand(),50))*(ux_set1(f,d)-lx_set1(f,d)));
  desp1{f,2}(d,r)=ly_set1(f,d)+(abs(weierstrass_function(rand(),50))*(uy_set1(f,d)-ly_set1(f,d)));
  desp1{f,3}(d,r)=lz_set1(f,d)+(abs(weierstrass_function(rand(),50))*(uz_set1(f,d)-lz_set1(f,d)));
  end
  end
end
% desp=dist_discrete(desp,desxy,dimorder,dima);% distance initialition accomplishment
% desp1=dist_discrete(desp1,desxy1,dimorder,dima);
%% change acceleration to distance +++++++++++++++++++++++ (multitasking)
for i=1:N,
     desp_tempx=[];desp_tempy=[];desp_tempz=[];
    for j=1:dimorder,
        desp_tempx=[desp_tempx desp{i,1}(j,:)];
        desp_tempy=[desp_tempy desp{i,2}(j,:)];
        desp_tempz=[desp_tempz desp{i,3}(j,:)];
    end
desp{i,1}=desp_tempx;desp{i,2}=desp_tempy;desp{i,3}=desp_tempz;
Salps_X(1,i)=desp(i,1);
Salps_X(2,i)=desp(i,2);
Salps_X(3,i)=desp(i,3);
end

for i=1:N,
     desp_tempx1=[];desp_tempy1=[];desp_tempz1=[];
    for j=1:dimorder,
        desp_tempx1=[desp_tempx1 desp1{i,1}(j,:)];
        desp_tempy1=[desp_tempy1 desp1{i,2}(j,:)];
        desp_tempz1=[desp_tempz1 desp1{i,3}(j,:)];
    end
desp1{i,1}=desp_tempx1;desp1{i,2}=desp_tempy1;desp1{i,3}=desp_tempz1;
Salps_X2(1,i)=desp1(i,1);
Salps_X2(2,i)=desp1(i,2);
Salps_X2(3,i)=desp1(i,3);
end

Salps_X=pos_bound_adjustment(Salps_X,lx_set,ux_set,ly_set,uy_set,lz_set,uz_set,N,dimorder,dimv);
Salps_X2=pos_bound_adjustment(Salps_X2,lx_set1,ux_set1,ly_set1,uy_set1,lz_set1,uz_set1,N,dimorder,dimv);
%% iteration begin
temp_number=0;
Archive_size=zeros(1,max_iter+1);
Archive_size(1)=0;
for iter=1:max_iter,
%    best_Solves1{k,iter}=Salps_X;
%    fitness_set1{k,iter}=Salps_fitness;
    c1 = 2*exp(-(4*iter/max_iter)^2); % Eq. (3.2) in the paper
    %% update again
f_1=@(xtt) xtt.^betaa.*exp(-xtt);f_1e=@(xtt) xtt.^betaa1.*exp(-xtt);f_1l=@(xtt) xtt.^betaa2.*exp(-xtt);
f_2=@(xtt) xtt.^(round((1+betaa)/2)-1).*exp(-xtt);f_2e=@(xtt) xtt.^(round((1+betaa1)/2)-1).*exp(-xtt);f_2l=@(xtt) xtt.^(round((1+betaa2)/2)-1).*exp(-xtt);
I_res1=integral(f_1,0,inf);I_res1e=integral(f_1e,0,inf);I_res1l=integral(f_1l,0,inf);
I_res2=integral(f_2,0,inf);I_res2e=integral(f_2e,0,inf);I_res2l=integral(f_2l,0,inf);
del_u=(I_res1*sin(pi*betaa/2)/(I_res2*betaa*2^((betaa-1)/2)))^(1/betaa);
del_ue=(I_res1e*sin(pi*betaa1/2)/(I_res2e*betaa1*2^((betaa1-1)/2)))^(1/betaa1);
del_ul=(I_res1l*sin(pi*betaa2/2)/(I_res2l*betaa2*2^((betaa2-1)/2)))^(1/betaa2);
del_v=1;
R_Gaussian_u=abs(del_u)*randn(1,dim-1);R_Gaussian_ue=abs(del_ue)*randn(1,dim-1);R_Gaussian_ul=abs(del_ul)*randn(1,dim-1);
R_Gaussian_v=abs(del_v)*randn(1,dim-1);
R_Gaussian_v=abs(R_Gaussian_v);
sss=R_Gaussian_u./R_Gaussian_v.^(1/betaa);sss1=R_Gaussian_ue./R_Gaussian_v.^(1/betaa1);sss2=R_Gaussian_ul./R_Gaussian_v.^(1/betaa2);
sss_record{iter}=sss;  
% c2=sin(pi)
   %% improved_section6
   wcc=rand()/2;wcc1=rand()/2;wcc2=rand()/2;
    c2=c2_lib(iter);
    %c2=rand(1,1);
    %c3=rand(1,1);
    c3=c3_lib(iter);
    c3e=c3_libe(iter);
    c3l=c3_libl(iter);
    c2m=c2_libm(iter);
    c3m=c3_libm(iter);
    c3me=c3_libe(iter);
    c3ml=c3_libl(iter);
    c2me=c2_libm(iter);
    c2ml=c2_libm(iter);
    for i=1:N, %Calculate all the objective values first
       % Salps_fitness(i,:)=ObjectiveFunction(Salps_X(:,i)');
        Salps_fitness(i,:)=ObjectiveFunction(Salps_X(:,i)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);%% function 
        Salps_fitness2(i,:)=ObjectiveFunction2(Salps_X2(:,i)',dim_length,penalty_constant*k,pos_xyz,acceler_bound); % multitasking
        if dominates(Salps_fitness(i,:),Food_fitness)
            Food_fitness=Salps_fitness(i,:);
            Food_position=Salps_X(:,i); % cell
        end
        % multitasking
        if dominates(Salps_fitness2(i,:),Food_fitness2)
            Food_fitnes2=Salps_fitness2(i,:);
            Food_position2=Salps_X2(:,i); 
        end
    end
  %% Problem Archive_X????
    [Archive_X, Archive_F, Archive_member_no]=UpdateArchive1(Archive_X, Archive_F, Salps_X, Salps_fitness, Archive_member_no);
    [Archive_X1, Archive_F1, Archive_member_no1]=UpdateArchive1(Archive_X1, Archive_F1, Salps_X2, Salps_fitness2, Archive_member_no1);% multitasking
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
        [Archive_X, Archive_F, Archive_mem_ranks, Archive_member_no]=HandleFullArchive(Archive_X, Archive_F, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
    end

    %% multitasking 
    if Archive_member_no1>ArchiveMaxSize
        Archive_mem_ranks1=RankingProcess(Archive_F1, ArchiveMaxSize, obj_no);
        [Archive_X1, Archive_F1, Archive_mem_ranks1, Archive_member_no1]=HandleFullArchive(Archive_X1, Archive_F1, Archive_member_no1, Archive_mem_ranks1, ArchiveMaxSize);
    else
        Archive_mem_ranks1=RankingProcess(Archive_F1, ArchiveMaxSize, obj_no);
    end
    Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
    Archive_mem_ranks1=RankingProcess(Archive_F1, ArchiveMaxSize, obj_no);% multitasking
    index=RouletteWheelSelection(1./Archive_mem_ranks);
    index1=RouletteWheelSelection(1./Archive_mem_ranks1);% multitasking
    if index==-1
        index=1;
    end
    %% multitasking
    if index1==-1
        index1=1;
    end
    Food_fitness=Archive_F(index,:);
    Food_position=Archive_X(index,:)';
    %% multitasking
    Food_fitness2=Archive_F1(index1,:);
    Food_position2=Archive_X1(index1,:)';
    %% 
    for g=1:size(Archive_X,1),
    if sum(Archive_X{g,dim})~=sum(1:1:dimorder+1)
    quit();
    else
    continue;
    end
    end
    %% recording population and fitness value 
    best_Solves1{k,iter}=Salps_X;
    fitness_set1{k,iter}=Salps_fitness;
    %% multitasking
    best_Solves2{k,iter}=Salps_X2;
    fitness_set2{k,iter}=Salps_fitness2;
    %% 
    %% multitasking update
    pop1_o=Salps_X;%CMOP_orignal
    pop2_o=Salps_X2;%MOP_orignal
 


if iter<=max_iter*iter_ratio
[Salps_X, final_disobj]=mtmssa_update_motion2_update2(fitness_set1,dimv,sss,best_Solves1,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2,c3,dimorder,ub_ud,lb_ud,dim_length,dim_total,wcc,iter,Archive_size,max_iter);   
temp_number=temp_number+1;
else
[Salps_X, final_disobj]=mssa_update_motion2(dim_total,dimv,dim,Salps_X,Food_position,c1,N,luratio,dima,pos_xyz,c2m,c3m,dimorder,ub_ud,lb_ud,dim_length);
end

early_number=0;later_number=0;
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       %% early phase-MOP
     if FES<iter_ratio*MAXFES*betta11 && iter<=iter_ratio*max_iter
         early_number=early_number+1;
       e_temp1=randperm(N,N/2);
       e_temp2=randperm(N,N/2);
       e_pop1=pop1_o(:,e_temp1);%CMOP_random
       e_pop2=pop2_o(:,e_temp2);%MOP_random
       e_pop1_u=Salps_X(:,e_temp1);%CMOP_updated
       Salps_X2=EMCMO_EARLY1(dima,c2me,c3me,max_iter,iter_ratio,luratio,iter,N,Food_position2,Salps_X2,dim,sss1,c1,c2,c3e,wcc1,dim_length,ub_ud1,lb_ud1,ub_total1,lb_total1,best_Solves2,fitness_set2,pos_xyz,dimv,dim_total,lx_set1,ly_set1,ux_set1,uy_set1,desxy1,dimorder,Archive_size);
       e_pop2_u=Salps_X2(:,e_temp2);%MOP_updated
       e_Tpop1=[pop1_o';e_pop1_u';e_pop2_u']';
       e_Tpop2=[pop2_o';e_pop2_u';e_pop1_u']';
       for f=1:size(e_Tpop1,2),
       e_Salps_Tpop1f(f,:)=ObjectiveFunction(e_Tpop1(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       e_Salps_Tpop2f(f,:)=ObjectiveFunction2(e_Tpop2(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       end
       %% environmental selection-Tpop1
       % [envir_con1,envir_con2]=sort(e_Salps_Tpop1f(:,1));[envir_con3,envir_con4]=sort(e_Salps_Tpop1f(:,2));
       % e_Salps_Tpop1f_sum=e_Salps_Tpop1f(:,1)+e_Salps_Tpop1f(:,2);% sum upe_Salps_Tpop1f
       %envir1=new_experimental_selection(e_Salps_Tpop1f,60);
       e_Salps_Tpop1f_sum=sqrt(e_Salps_Tpop1f(:,1).^2+e_Salps_Tpop1f(:,2).^2);
       [envir_con_sum1, envir_con_sum2]=sort(e_Salps_Tpop1f_sum);% sort 
         Salps_X=e_Tpop1(:,envir_con_sum2(1:N));
       %% environmental selection-Tpop2
       % [envir_con5,envir_con6]=sort(e_Salps_Tpop2f(:,1));[envir_con7,envir_con8]=sort(e_Salps_Tpop2f(:,2));
       % e_Salps_Tpop2f_sum=e_Salps_Tpop2f(:,1)+e_Salps_Tpop2f(:,2);% sum upe_Salps_Tpop2f
       e_Salps_Tpop2f_sum=sqrt(e_Salps_Tpop2f(:,1).^2+e_Salps_Tpop2f(:,2).^2);
      % envir2=new_experimental_selection(e_Salps_Tpop1f,60);
       [envir_con_sum3, envir_con_sum4]=sort(e_Salps_Tpop2f_sum);% sort  
       Salps_X2=e_Tpop2(:,envir_con_sum4(1:N));
     elseif FES>=iter_ratio*MAXFES*betta11 && iter<=iter_ratio*max_iter
     %% Later phase-MOP
      later_number=later_number+1;
       for f=1:size(pop1_o,2),
       l_Salps_pop1o(f,:)=ObjectiveFunction(pop1_o(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       l_Salps_pop2o(f,:)=ObjectiveFunction2(pop2_o(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       end
       l_temp1=Tournament_Selection(l_Salps_pop1o,N/2);
       l_temp2=Tournament_Selection(l_Salps_pop2o,N/2);
       if ~(isnumeric(l_temp1) && isreal(l_temp1) && all(l_temp1 > 0) && all(mod(l_temp1, 1) == 0)) && ~islogical(l_temp1)
       %error('SafeArrayAccess:InvalidIndex', 'Array index must be a positive integer or logical value.');
       error('SafeArrayAccess:InvalidIndex', ...
        ['Array index must be a positive integer or logical value. Value provided: ', mat2str(l_temp1),num2str(length(l_temp1))]);
       end
       if ~(isnumeric(l_temp2) && isreal(l_temp2) && all(l_temp2 > 0) && all(mod(l_temp2, 1) == 0)) && ~islogical(l_temp2)
       %error('SafeArrayAccess:InvalidIndex', 'Array index must be a positive integer or logical value.');
           error('SafeArrayAccess:InvalidIndex', ...
          ['Array index must be a positive integer or logical value. Value provided: ', mat2str(l_temp2),num2str(length(l_temp2))]);
       end
       l_pop1=pop1_o(:,l_temp1);
       l_pop2=pop2_o(:,l_temp2);
       l_pop1_u=Salps_X(:,l_temp1);
      %[lx_set1,ux_set1,ly_se1t,uy_set1]=upper_lower_update(pos_xyz,Salps_X2,dima); 
       Salps_X2=EMCMO_LATER1(dima,c2ml,c3ml,max_iter,iter_ratio,luratio,iter,N,Food_position2,Salps_X2,dim,sss2,c1,c2,c3l,wcc2,dim_length,ub_ud1,lb_ud1,ub_total1,lb_total1,best_Solves2,fitness_set2,pos_xyz,dimv,dim_total,lx_set1,ly_set1,ux_set1,uy_set1,desxy1,dimorder,Archive_size);
       l_pop2_u=Salps_X2(:,l_temp2);
       l_Tpop1=[pop1_o';l_pop1_u']';% 3/2 N
       l_Tpop2=[pop2_o';l_pop2_u']';
       for f=1:size(l_Tpop1,2),
       l_Salps_Tpop1o(f,:)=ObjectiveFunction2(l_Tpop1(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       l_Salps_Tpop2o(f,:)=ObjectiveFunction(l_Tpop2(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       end
       %........................................Trp1 and Trp2
       [l_Tp1,l_Tp2]=ENVIRONMENT_SELECTION_JUDGE(l_Salps_Tpop1o,l_Salps_Tpop2o,l_Tpop1,l_Tpop2,N,size(pop1_o,2),size(l_pop1_u,2));
       for f=1:size(l_Tp1,2),
       l_Salps_Tp1(f,:)=ObjectiveFunction(l_Tp1(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       l_Salps_Tp2(f,:)=ObjectiveFunction2(l_Tp2(:,f)',dim_length,penalty_constant*k,pos_xyz,acceler_bound);
       end
       [Salps_X,Salps_X2]=ENVIRONMENT_SELECTION(l_Salps_Tp1,l_Salps_Tp2,l_Tp1,l_Tp2,N);
     else
         Salps_X=Salps_X;
     end
     iter_early_later(iter,:)=[early_number,later_number];

     FES=FES+N;


[lx_set,ux_set,ly_set,uy_set,lz_set,uz_set,desxy]=upper_lower_update1(pos_xyz,Salps_X,dim,N,luratio); % update upper and lower bounds
for i=1:N
         for v=4:4,
            for l=1:length(final_disobj{i}),
                if final_disobj{i}(l)==1
        Flag4ub=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)>ub_ud(1,v);
        Flag4lb=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)<lb_ud(1,v);
        Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)=(Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l).*(~(Flag4ub+Flag4lb)))+ub_ud(1,v).*Flag4ub+lb_ud(1,v).*Flag4lb;  
                else
        Flag4ub=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)>ub_ud(2,v);
        Flag4lb=Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)<lb_ud(2,v);
        Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l)=(Salps_X{v,i}(1+(dim_length(v)*(l-1)):dim_length(v)*l).*(~(Flag4ub+Flag4lb)))+ub_ud(2,v).*Flag4ub+lb_ud(2,v).*Flag4lb;     
                end
             end
         end
        for v=1:3,
            if v==1
                for l=1:dimorder,
        Flag4ub=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))>max(ux_set(i,l),lx_set(i,l));
        Flag4lb=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))<min(lx_set(i,l),ux_set(i,l));
        Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*(~(Flag4ub+Flag4lb)))+max(ux_set(i,l),lx_set(i,l)).*Flag4ub+min(lx_set(i,l),ux_set(i,l)).*Flag4lb;  
        %Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*desxy{i,1}(l)./sum(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)));
                end
            elseif v==2
                for l=1:dimorder,
        Flag4ub=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))>max(uy_set(i,l),ly_set(i,l));
        Flag4lb=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))<min(ly_set(i,l),uy_set(i,l));
        Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*(~(Flag4ub+Flag4lb)))+max(uy_set(i,l),ly_set(i,l)).*Flag4ub+min(ly_set(i,l),uy_set(i,l)).*Flag4lb;  
                %Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*desxy{i,2}(l)./sum(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)));
                end
            elseif v==3
                for l=1:dimorder,
        Flag4ub=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))>max(uz_set(i,l),lz_set(i,l));
        Flag4lb=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))<min(lz_set(i,l),uz_set(i,l));
        Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*(~(Flag4ub+Flag4lb)))+max(uz_set(i,l),lz_set(i,l)).*Flag4ub+min(lz_set(i,l),uz_set(i,l)).*Flag4lb;  
                %Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v))=Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)).*desxy{i,2}(l)./sum(Salps_X{v,i}(1+(l-1)*dim_length(v):l*dim_length(v)));
                end
            end
        end
end
    disp(['iteration:', num2str(iter), ', non-dominated solutions number:', num2str(Archive_member_no)]);
Salps_X=pos_bound_adjustment(Salps_X,lx_set,ux_set,ly_set,uy_set,lz_set,uz_set,N,dimorder,dimv);
Salps_X2=pos_bound_adjustment(Salps_X2,lx_set1,ux_set1,ly_set1,uy_set1,lz_set1,uz_set1,N,dimorder,dimv);
best_Solves{k,iter}=Archive_X;
%fitness_set{k,iter} = Archive_F;
fitness_set{k,iter}=Archive_F;
fitness_set1{k,iter}=Archive_F1;
salps_state{k,iter}=Salps_X;
sfalps_state{k,iter}=Salps_fitness;
Archive_size(iter+1)=size(Archive_F,1);
end
% best_Solves{k,i} = Archive_X;
% fitness_set{k,i} = Archive_F;
end
% file_name = ['C:\Users\Administrator\Desktop\MSSA_MSSA_update\mtmssa', num2str(cycle_number+10), '.mat'];
% eval(['save(''', file_name, ''', ''fitness_set'',''best_solvers'');']);
end
% save('C:\Users\Administrator\Desktop\MSSA_MSSA_update\mtmssa12.mat','fitness_set','best_Solves');
% figure
% Draw_ZDT1();
% hold on
% %plot(Archive_F(:,1),Archive_F(:,2),'ro','MarkerSize',8,'markerfacecolor','k');
% figure
% plot(fitness_set{1,500}(1:end,1),fitness_set{1,500}(1:end,2),'ro','MarkerSize',8,'markerfacecolor','k');
% %legend('True PF','Obtained PF');
% title('MSSA');