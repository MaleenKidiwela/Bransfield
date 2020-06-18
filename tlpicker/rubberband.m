function [xmin, xmax, ymin, ymax] = rubberband(handleFigure)

if nargin>=1
  figure(handleFigure)
end

waitforbuttonpress;
point1 = get(gca,'CurrentPoint');  
finalRect = rbbox;                   
point2 = get(gca,'CurrentPoint');    
point1 = point1(1,1:2);              
point2 = point2(1,1:2);  

xmin = min(point1(1),point2(1));
xmax = max(point1(1),point2(1));
ymin = min(point1(2),point2(2));
ymax = max(point1(2),point2(2));  