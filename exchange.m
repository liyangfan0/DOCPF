function box2 = exchange( box1 )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
x=box1(1);
y=box1(2);
w=box1(3);
h=box1(4);
box2=[x-0.5*w,y-0.5*h,x+0.5*w,y+0.5*h];

end

