function output=APCE(response)
%计算APCE的值
%response:响应值，是一个二维高斯分布
%output:输出一个值
ymax=max(response(:));
ymin=min(response(:));
a=(ymax-ymin)^2;
b=(response-ymin).^2;
output=a/sum(b(:));

