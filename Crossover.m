function [yy1,yy2]=Crossover(x1,x2)
    for i=1:size(x1,1)-1,
        % x1=pop(1).Position;
        xx1=x1{i};xx2=x2{i};
    %alpha=rand(size(xx1)); 
    % gamma=0.5;
    % alpha=-gamma+(1+2*gamma)*rand(1,length(xx1));
    alpha=rand(size(xx1));
    for j=1:length(alpha),
    if  alpha(j)<=0.5
        B(j)=(2*alpha(j))^(1/(sqrt(2*pi)+1));
    else
        B(j)=(2*(1-alpha(j)))^(1/(sqrt(2*pi)+1)); 
    end
    end
    %alpha=alpha/1000;
     %alpha=rand(size(xx1)); 
    y1=((1-B).*xx1+(1+B).*xx2)/2;
    y2=((1-B).*xx2+(1+B).*xx1)/2;
    % y1=alpha.*xx1+(1-alpha).*xx2;
    % y2=alpha.*xx2+(1-alpha).*xx1;
    x1{i}=y1;x2{i}=y2;
    end
    % r=rand();
  % if r<delta1
  % x1{5}=randperm(length(x1{5}));
  % x2{5}=randperm(length(x2{5}));
  % elseif r>=delta1 && r<delta2
  %   x1{5}=x1{5};x2{5}=x2{5};
  % else
  %     x1{5}=best_position{5};
  %     x2{5}=best_position{5};
  % end
  yy1=x1;yy2=x2;
end