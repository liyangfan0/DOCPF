function IoU = compute_IoU(box1, box2)
%COMPUTE_IOU Is compute the two region overlap area.
%    
%   ************************
%   *                      *
%   *      (x_a,y_a)******************
%   *          *           *         *
%   *          *           *         *
%   *          *           *         *
%   *******************(x_b,y_b)     *
%              *                     *
%              *                     *
%              ***********************
% ‰»Îcx,cy,w,h
box1=exchange(box1);
box2=exchange(box2);
width1=abs(box1(3)-box1(1));
height1=abs(box1(2)-box1(4));
width2=abs(box2(3)-box2(1));
height2=abs(box2(2)-box2(4));
xmax=max([box1(1),box1(3),box2(1),box2(3)]);
ymax=max([box1(2),box1(4),box2(2),box2(4)]);
xmin=min([box1(1),box1(3),box2(1),box2(3)]);
ymin=min([box1(2),box1(4),box2(2),box2(4)]);
W=xmin+width1+width2-xmax;
H=ymin+height1+height2-ymax;
if W<=0 || H<=0
    IoU=0;
else
    iou_area=W*H;
    box1_area=width1*height1;
    box2_area=width2*height2;
    IoU=iou_area/(box1_area+box2_area-iou_area);
end

end




