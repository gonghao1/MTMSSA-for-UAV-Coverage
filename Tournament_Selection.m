function sell1=Tournament_Selection(Fit,N),
% sell1=[];
% Sel1=[];
% for i = 1:size(Fit,1)
%     temp = randi([1 size(Fit,1)],1,2);   
%     if temp(1)==0 
%         temp(1)=1;
%     elseif temp(2)==0
%         temp(2)=1;
%     end
%     if Fit(temp(1),1)<=Fit(temp(2),1) && Fit(temp(1),2)<=Fit(temp(2),2)
%         Sel1(i)=temp(1);
%     elseif Fit(temp(1),1)<=Fit(temp(2),1) && Fit(temp(1),2)>Fit(temp(2),2)
%         Sel1(i)=temp(1);
%     elseif Fit(temp(1),1)>Fit(temp(2),1) && Fit(temp(1),2)>Fit(temp(2),2)
%         Sel1(i)=temp(2);
%     elseif Fit(temp(1),1)>Fit(temp(2),1) && Fit(temp(1),2)<=Fit(temp(2),2)
%         Sel1(i)=temp(2);
%     end
% end
% %Sel1=Sel1';
% while length(sell1)~=N,
%    sell1=[sell1 mode(Sel1)];
%    %a=find(Sel1==sell1);
%    for i=1:length(sell1),
%    Sel1(find(Sel1==sell1(i)))=[];
%    end
% end
% if length(find(sell1==0))~=0
% sell1(find(sell1==0))=1;
% end
% 检查输入
    if isempty(Fit) || all(all(isnan(Fit)))
        error('Fitness matrix is empty or all values are NaN.');
    end

    % 初始化选择数组
    Sel1 = zeros(size(Fit, 1), 1); % 用0初始化
    for i = 1:size(Fit, 1)
        temp = randi([1 size(Fit, 1)], 1, 2); % 随机选两个个体
        
        % 锦标赛选择逻辑
        if isnan(Fit(temp(1), 1)) || isnan(Fit(temp(1), 2)) || isnan(Fit(temp(2), 1)) || isnan(Fit(temp(2), 2))
            % 如果选中的个体的适应度是NaN，选择一个未在Sel1出现的最小正整数
            Sel1(i) = find(~ismember(1:max(Sel1)+1, Sel1), 1, 'first');
        elseif (Fit(temp(1), 1) <= Fit(temp(2), 1) && Fit(temp(1), 2) <= Fit(temp(2), 2)) || ...
               (Fit(temp(1), 1) <= Fit(temp(2), 1) && Fit(temp(1), 2) > Fit(temp(2), 2))
            Sel1(i) = temp(1);
        else
            Sel1(i) = temp(2);
        end
    end
    
    % 清除初始化时的0
    Sel1 = Sel1(Sel1 > 0);

    % 获取唯一的选择结果
    unique_sel1 = unique(Sel1);
    
    % 确保有足够的唯一元素进行选择
    if length(unique_sel1) < N/2
        error('Not enough unique selected elements to meet the requirement.');
    end

    % 从唯一值中随机选择N/2个元素
    sell1 = randsample(unique_sel1, N/2);

% final_fit{k}=sell1;
end
% trr=[];
% for i=1:length(final_fit),
%     if length(find(final_fit{i}==0))==0
%     trr(i)=1;
%     else
%         trr(i)=0;
%     end
% end