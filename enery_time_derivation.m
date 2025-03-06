function record_final=enery_time_derivation(c,pos_final),
ObjectiveFunction=@ZDT5_rr;
% acce_bound=[-3,3,-1,1];
%c=fitness4_mt2;
for i=1:length(c),
xtee=c{i};
dsize=length(xtee);
pos_tt=pos_final{i};
%xteee=xtee{dsize};
M=length(xtee{1,1})/(length(xtee{1,5})-1);
M=[M,M,M,M];
er=[];tr=[];fr=[];
for k=1:size(xtee,1),
er_tr=ObjectiveFunction(xtee,M,pos_tt);
% er(k)=er_tr(1);
% tr(k)=er_tr(2);
fr=[fr;er_tr];
end
%fr=[er;tr]';
record_final{i}=fr;
end
end