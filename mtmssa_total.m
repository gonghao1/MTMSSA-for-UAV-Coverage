function [best_Solves,fitness_set,fitness_set1]=mtmssa_total(x1_fin1,y1_fin1,UAV_number,UAV_onboard,max_iter,N,ty,pos_final),
ObjectiveFunction=@ZDT3; % 目标函数
dim=2;dimv=501;dima=dimv-1;
%dimorder=10;% ？？？变量维度
lb=0;% 变量下界
ub=1;% 变量上界
obj_no=2;%目标函数数量
pos_final=pos_final(ty);
region_number=length(y1_fin1);% 可以自行调整
for i=1:region_number,
fitness_tem=x1_fin1{i};
fitness_tem1=y1_fin1{i};
dft_1=sqrt(fitness_tem(:,1).^2+fitness_tem(:,2).^2);
fitness_num=find(dft_1==min(dft_1));
fitness_num=fitness_num(1);
%find(fitness_tem(:,1)==min(fitness_tem(:,1)));
fitness_tem1_opt=fitness_tem1(fitness_num,:);
fitness_tem1_opt_f{i}=fitness_tem1_opt;
xx=fitness_tem1_opt{end};
reg_num_sin(i)=length(fitness_tem1_opt{end});
fitness_num_record(i)=fitness_num;
fitness_path_record{i}=fitness_tem1_opt;
end
%save('202407311.mat',"fitness_path_record");
pos_xyz=[];
for i=1:region_number,
pos_xyz=[pos_xyz;pos_final{i}];
end
new_time_energy=enery_time_derivation(fitness_tem1_opt_f,pos_final);
fnte_e=[];fnte_t=[];
for i=1:length(new_time_energy),
nte=new_time_energy{i};
nte_e=nte{1};nte_t=nte{2};
sepe_rec{i}=nte_e;sept_rec{i}=nte_t;
fnte_e=[fnte_e nte_e];
fnte_t=[fnte_t nte_e];
end

region_point_number=size(pos_xyz,1);
centerp=mean(pos_xyz(:,1:2));
posx_in=centerp(1);posy_in=centerp(2);posz_in=55;
pos_center=[posx_in,posy_in,posz_in];
cp_e=[];cp_t=[];pc_e=[];pc_t=[];
% for i=1:region_point_number,
% array_cp=[pos_center;
%     pos_xyz(i,:)];
% array_pc=[pos_xyz(i,:)
%     pos_center];
% [a_cp,b_cp]=mtmssa_test_twopoint(0.6,0,5,array_cp,1);
% [a_pc,b_pc]=mtmssa_test_twopoint(0.6,0,5,array_cp,1);
% optimal_cp_num=find(a_cp{end}(:,1)==min(a_cp{end}(:,1)));
% optimal_cp=b_cp{end}(optimal_cp_num,:);
% optimal_pc_num=find(a_pc{end}(:,1)==min(a_pc{end}(:,1)));
% optimal_pc=b_pc{end}(optimal_pc_num,:);
% cp_et=enery_time_derivation1(optimal_cp,array_cp);
% pc_et=enery_time_derivation1(optimal_pc,array_pc);
% cp_e=[cp_e cp_et{1}];cp_t=[cp_t cp_et{2}];
% pc_e=[pc_e cp_et{1}];pc_t=[pc_t pc_et{2}];
% fcp_a{i}=a_cp;fcp_b{i}=b_cp;
% fpc_a{i}=a_pc;fpc_b{i}=b_pc;
% end

%save('20240911.mat',"cp_t","cp_e","pc_e","pc_t");

%load('202409103.mat');

load('2024091231.mat');

ArchiveMaxSize=800;
epslison=1;
max_runtimes=1;
penalty_constant=2;

%% algorithm test !!!!!!
%segments_UAV=fixedPartition(desxyz, UAV_onboard, UAV_number);% resure UAV number

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
%% algorithm_begin
for k=1:max_runtimes,
%% archive
%Archive_X=zeros(100,dim);
for i=1:dim,
for j=1:ArchiveMaxSize,
Archive_X{i,j}=zeros(1,UAV_number);
end
end
Archive_X=Archive_X';
Archive_F=ones(ArchiveMaxSize,obj_no)*inf;
Archive_F1=ones(ArchiveMaxSize,obj_no)*inf;
Archive_member_no=0;
%r=(ub-lb)/2;
%V_max=(ub(1)-lb(1))/10;
Food_fitness=inf*ones(1,obj_no);
Food_position=zeros(1,1);
%Salps_X=initialization(N,dim,ub,lb);
%% algorithm test !!!
for m=1:N,
Salps_X{2,m}=splitPositiveIntegerMIntoKParts(region_point_number, UAV_number);
Salps_X{1,m}=randperm(region_number);
end
% xx=Salps_X(:,1);
%%　improved_section7 levy flight mechanism
for iter=1:max_iter,
    c1 = 2*exp(-(4*iter/max_iter)^2); % Eq. (3.2) in the paper
    c11 = 2*exp(-(4*iter/max_iter)^2); % Eq. (3.2) in the paper
    % c1 = c2_lib(iter); % Eq. (3.2) in the paper
    % c11 =c2_lib(iter);% Eq. (3.2) in the paper
   % c2=sin(pi)
   %% improved_section6
   % wcc=rand()/2;
   %  c2=c2_lib(iter);
   %  c3=c3_lib(iter);
   for i=1:N, %Calculate all the objective values first
        otemp=ObjectiveFunction(Salps_X(:,i)',sepe_rec,sept_rec,pc_e,pc_t,cp_e,cp_t,UAV_onboard,penalty_constant*k,reg_num_sin);%% function 
        Salps_fitness(i,:)=otemp{1};
        Salps_fitness1(i,:)=otemp{2};
        if dominates(Salps_fitness(i,:),Food_fitness)
            Food_fitness=Salps_fitness(i,:);
            Food_position=Salps_X{i}; % cell
        end
   end
  %% Problem Archive_X????
    [Archive_X, Archive_F,Archive_F1, Archive_member_no]=UpdateArchive_total(Archive_X, Archive_F, Archive_F1,Salps_X, Salps_fitness,Salps_fitness1,Archive_member_no);
    
    if Archive_member_no>ArchiveMaxSize
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
        [Archive_X, Archive_F, Archive_F1,Archive_mem_ranks, Archive_member_no]=HandleFullArchive_total(Archive_X, Archive_F, Archive_F1, Archive_member_no, Archive_mem_ranks, ArchiveMaxSize);
    else
        Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);
    end

    Archive_mem_ranks=RankingProcess(Archive_F, ArchiveMaxSize, obj_no);

    index=RouletteWheelSelection(1./Archive_mem_ranks);
    if index==-1
        index=1;
    end
    Food_fitness=Archive_F(index,:);
    Food_position=Archive_X(index,:)';
    %% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    % for g=1:size(Archive_X,1),
    % if sum(Archive_X{g,5})~=sum(1:1:dimorder+1)
    % quit();
    % else
    % continue;
    % end
    % end
    best_Solves1{k,iter}=Salps_X;
    fitness_set1{k,iter}=Salps_fitness;
    fitness_set2{k,iter}=Salps_fitness1;
    %% 
 %% Update salps 
    for i=1:N,
             c1=rand;
            % f_order1=PMX_function_or(Food_position{1},Salps_X{1,i},c1);
            f_order1=randperm(length(Food_position{1}));
             Salps_X{1,i}=f_order1;
             c11=rand();
             f_order2=pmx_mut_function(Food_position{2},Salps_X{2,i},c11);
             Salps_X{2,i}=f_order2;
    end
     
    disp(['At the iteration ', num2str(iter), ' there are ', num2str(Archive_member_no), ' non-dominated solutions in the archive']);
best_Solves{k,iter}=Archive_X;
fitness_set{k,iter}=Archive_F;
fitness_set1{k,iter}=Archive_F1;
end
end