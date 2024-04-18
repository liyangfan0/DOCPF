function flag=flag_count(S1,S2)
%S1是记录的APCE值，S2是记录的Fmax的值
S{1}=[mean(S2),mean(S1)]; %当前APCE值和Fmax的值与历史平均值相距甚远，此时判断为遮挡
S{2}=[mean(S2(1:end-1)),mean(S1(1:end-1))]; %有时候旋转也会导致响应值不高，此时为了区分好旋转问题，需进一步处理
%S{3}=[mean(S2(1:end-2)),mean(S1(1:end-2))];

%根据置信度求flag的值%0.35,0.35,0.1
if (S2(end)/S{2}(1)<0.35) &&(S1(end)/S{2}(2)<=0.35) && S2(end)>=0.1 %旋转的情况
    flag=1;
else
    if(S2(end)/S{2}(1)>=0.35) && (S1(end)/S{2}(2)>=0.35)%&&(S2(end-1)/S{2}(1)>=0.3)
        flag=1;
    else
        flag=2;
    end
end
%else
%   if (S2(end)/S{1}(1)<0.35) &&(S1(end)/S{1}(2)<=0.35) && S2(end)>=0.1 %旋转的情况
%       flag=1;
%   else
%       if(S2(end)/S{1}(1)>=0.35) &&(S1(end)/S{1}(2)>=0.35)&&(S2(end-1)/S{2}(1)>=0.3)
%           flag=1;
%       else
%           flag=2;
%       end
%   end
%end