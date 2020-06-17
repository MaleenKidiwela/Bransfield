p1 = load('T5_22_01_19_EDIT.asvp');
p3 = load('T5_27_01_19.asvp');
p4 = load('T5_30_01_19.asvp');


Average_OBS_Depth = (1455+1496+1345+1336+1460+1239+1087+1370+1092+1025+1367+1157+916+1433+1475)/15;
%lower limit is at 15 m (source depth) and upper limit is at 1283m(average
%obs depth)

h = round(Average_OBS_Depth)-15;
max_D =round(Average_OBS_Depth);
min_D = -1*srEvent.elevation(1)*1000;

%interpolating in to a grid of 1m depth intervals
x = p1(:,1);
y = p1(:,2);
x1 = linspace(min_D,max_D,h);
y1 = interp1(x,y,x1,'linear');
clear x y
x = p3(:,1);
y = p3(:,2);
x2 = linspace(min_D,max_D,h);
y2 = interp1(x,y,x2,'linear');
clear x y
x = p4(:,1);
y = p4(:,2);
x3 = linspace(min_D,max_D,h);
y3 = interp1(x,y,x3,'linear');
clear x y


mean1 = sum(y1)/h;
mean2= sum(y2)/h;
mean3  = sum(y3)/h;

v=(mean1+mean2+mean3)/3;

[XVT5_WOA,z1]=getlev(-62.4167,-58.3833,'c');
[XVT6_WOA,z2]=getlev(-62.4083,-58.5333,'c');
[XVT2_WOA,z3]=getlev(-62.4700,-58.2213,'c');


xx = z1;
yy = XVT5_WOA;
xx1 = linspace(min_D,max_D,h);
yy1 = interp1(xx,yy,xx1,'linear');
clear x y
xx = z2;
yy = XVT6_WOA;
xx2 = linspace(min_D,max_D,h);
yy2 = interp1(xx,yy,xx2,'linear');
clear x y
xx = z3;
yy = XVT2_WOA;
xx3 = linspace(min_D,max_D,h);
yy3 = interp1(xx,yy,xx3,'linear');

mean11 = sum(yy1)/h;
mean22= sum(yy2)/h;
mean33  = sum(yy3)/h;

vv=(mean11+mean22+mean33)/3;
% 
% figure(80)
% plot(y1,-1*x1,'LineWidth',2)
% hold on
% plot(y2,-1*x2,'LineWidth',2)
% plot(y3,-1*x3,'LineWidth',2)
% plot(XVT5_WOA(1:22),-1*z1(1:22),':','LineWidth',2)% plotting depths upto 1300
% plot(XVT6_WOA(1:22),-1*z2(1:22),':','LineWidth',2)% plotting depths upto 1300
% plot(XVT2_WOA(1:22),-1*z3(1:22),':','LineWidth',2)% plotting depths upto 1300
% xline((mean1+mean2+mean3)/3,'r')
% xline(vv(1),'b')
% legend('XBT2','XBT5','XBT6','XBT5_W_O_A','XBT6_W_O_A','XBT2_W_O_A','average_X_B_T','average_W_O_A')
% 
% 
% 

clear Average_OBS_Depth mean1 mean2 mean3 p1 p3 p4 x1 x2 x3 y1 y2 y3

