function [a,b]=ENVIRONMENT_SELECTION(c,d,e,f,N),
       %c=l_Salps_Tp1;d=l_Salps_Tp2;e=l_Tp1;f=l_Tp2;
       %[envir_con1,envir_con2]=sort(c(:,1));[envir_con3,envir_con4]=sort(c(:,2));
       %e_Salps_Tpop1f_sum=c(:,1)+c(:,2);% sum upe_Salps_Tpop1f

       e_Salps_Tpop1f_sum=sqrt(c(:,1).^2+c(:,2).^2);% sum upe_Salps_Tpop1f
       %envir1=new_experimental_selection(c,60);
       [envir_con_sum1, envir_con_sum2]=sort(e_Salps_Tpop1f_sum);% sort 
       % envir_select_temp1=unique([envir_con4(1:N) envir_con2(1:N)]);
       % envir_con_sum2(envir_select_temp1)=[];% remove
      %  if length(envir_select_temp1)==N
      %  a=e(:,envir_select_temp1);
      %  elseif length(envir_select_temp1)>N
      %  a=e(:,envir_select_temp1(1:N));
      %  else
      %  envir_con2(envir_select_temp1)=[];
      % % envir_select_temp1=sort([envir_select_temp1 envir_con4(N+1:2*N-length(envir_select_temp1))]);
      %  envir_select_temp1=sort([envir_select_temp1 envir_con_sum2(1:N-length(envir_select_temp1))]);% replace 'envir_con6' with 'envir_con_sum2'  
      %  a=e(:,envir_select_temp1);
      %  end    
      a=e(:,envir_con_sum2(1:N));

       % [envir_con5,envir_con6]=sort(d(:,1));[envir_con7,envir_con8]=sort(d(:,2));
       % e_Salps_Tpop2f_sum=d(:,1)+d(:,2);% sum upe_Salps_Tpop2f
      
       e_Salps_Tpop2f_sum=sqrt(d(:,1).^2+d(:,2).^2);
       %envir2=new_experimental_selection(d,60);
       [envir_con_sum3, envir_con_sum4]=sort(e_Salps_Tpop2f_sum);% sort 
       % envir_select_temp2=unique([envir_con8(1:N) envir_con6(1:N)]);
       % envir_con_sum4(envir_select_temp2)=[];% remove 
       % if length(envir_select_temp2)==N
       % b=f(:,envir_select_temp2);
       % elseif length(envir_select_temp2)>N
       % b=f(:,envir_select_temp2(1:N));
       % else
       % envir_con6(envir_select_temp2)=[];
       % %envir_select_temp2=sort([envir_select_temp2 envir_con8(N+1:2*N-length(envir_select_temp2))]);
       % envir_select_temp2=sort([envir_select_temp2 envir_con_sum4(1:N-length(envir_select_temp2))]);% replace 'envir_con6' with 'envir_con_sum2'
       % b=f(:,envir_select_temp2);
       % end
       b=f(:,envir_con_sum4(1:N));
       %b=[f';trp1']';
end