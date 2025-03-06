function [Archive_X_updated, Archive_F_updated, Archive_member_no]=UpdateArchive1(Archive_X, Archive_F, Particles_X, Particles_F, Archive_member_no)
 Archive_X_updated={
     };
 Archive_F_updated=[];
 re_frepeat=[];
Archive_X_temp=[Archive_X; Particles_X'];% 维度由种群数量和archive  300x5
Archive_F_temp=[Archive_F; Particles_F];% 300x2
% isDominated = false(size(Archive_F_temp, 1), 1);
% for i = 1:size(Archive_F_temp, 1)
%     for j = i+1:size(Archive_F_temp, 1)
%         if dominates(Archive_F_temp(i, :), Archive_F_temp(j, :))
%             isDominated(j) = true;
%         elseif dominates(Archive_F_temp(j, :), Archive_F_temp(i, :))
%             isDominated(i) = true;
%             break; 
%         end
%     end
% end
findnan33=imag(Archive_F_temp(:,1));findnan34=imag(Archive_F_temp(:,2));
find_33=find(findnan33~=0);find_34=find(findnan34~=0);
find_3334=unique([find_33 find_34]);
Archive_X_temp(find_3334,:)=[];
Archive_F_temp(find_3334,:)=[];
% 初始化支配标志
isDominated = false(size(Archive_F_temp, 1), 1);
for i = 1:size(Archive_F_temp, 1)
    for j = 1:size(Archive_F_temp, 1)
        if i ~= j
            if dominates(Archive_F_temp(i, :), Archive_F_temp(j, :))
                isDominated(j) = true;
            elseif dominates(Archive_F_temp(j, :), Archive_F_temp(i, :))
                isDominated(i) = true;
                break; 
            end
        end
    end
    if isDominated(i) 
        continue;
    end
end
Archive_X_updated = Archive_X_temp(~isDominated, :);
Archive_F_updated = Archive_F_temp(~isDominated, :);
[findinf, ~]=find((Archive_F_updated==Inf) | (Archive_F_updated==-Inf));
findnan1=find(isnan(Archive_F_updated(:,1))==1);
findnan2=find(isnan(Archive_F_updated(:,2))==1);
findnan3=imag(Archive_F_updated(:,1));
findnan3=find(findnan3~=0);
findnan=unique([findnan1' findnan2']);
findinf=unique([findinf' findnan findnan3]);
if ~isempty(findinf)
Archive_X_updated(findinf,:)=[];
Archive_F_updated(findinf,:)=[];
end
% delate repeated items
[frepeat_1,frepeat_2]=unique(Archive_F_updated(:,1));
tempnum=0;
for i=1:length(frepeat_2),
if (imag(Archive_F_updated(i,1))*imag(Archive_F_updated(i,2)))==0
    tempnum=tempnum+1;
re_frepeat(tempnum)=frepeat_2(i);
end
end
Archive_X_updated = Archive_X_updated(re_frepeat,:);
Archive_F_updated = Archive_F_updated(re_frepeat,:);
Archive_member_no=size(Archive_F_updated,1);
end