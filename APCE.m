function output=APCE(response)
%����APCE��ֵ
%response:��Ӧֵ����һ����ά��˹�ֲ�
%output:���һ��ֵ
ymax=max(response(:));
ymin=min(response(:));
a=(ymax-ymin)^2;
b=(response-ymin).^2;
output=a/sum(b(:));

