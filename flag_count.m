function flag=flag_count(S1,S2)
%S1�Ǽ�¼��APCEֵ��S2�Ǽ�¼��Fmax��ֵ
S{1}=[mean(S2),mean(S1)]; %��ǰAPCEֵ��Fmax��ֵ����ʷƽ��ֵ�����Զ����ʱ�ж�Ϊ�ڵ�
S{2}=[mean(S2(1:end-1)),mean(S1(1:end-1))]; %��ʱ����תҲ�ᵼ����Ӧֵ���ߣ���ʱΪ�����ֺ���ת���⣬���һ������
%S{3}=[mean(S2(1:end-2)),mean(S1(1:end-2))];

%�������Ŷ���flag��ֵ%0.35,0.35,0.1
if (S2(end)/S{2}(1)<0.35) &&(S1(end)/S{2}(2)<=0.35) && S2(end)>=0.1 %��ת�����
    flag=1;
else
    if(S2(end)/S{2}(1)>=0.35) && (S1(end)/S{2}(2)>=0.35)%&&(S2(end-1)/S{2}(1)>=0.3)
        flag=1;
    else
        flag=2;
    end
end
%else
%   if (S2(end)/S{1}(1)<0.35) &&(S1(end)/S{1}(2)<=0.35) && S2(end)>=0.1 %��ת�����
%       flag=1;
%   else
%       if(S2(end)/S{1}(1)>=0.35) &&(S1(end)/S{1}(2)>=0.35)&&(S2(end-1)/S{2}(1)>=0.3)
%           flag=1;
%       else
%           flag=2;
%       end
%   end
%end