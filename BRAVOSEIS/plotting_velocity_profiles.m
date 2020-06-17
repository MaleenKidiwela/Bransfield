% here we are interested in generating velocity profiles and identifying a
% good constant velocity 

% velocity sounding for the area of interest is from the files
% for XVT2, XVT4, XVT5, XVT6
%
p1 = load('T5_22_01_19_EDIT.asvp');
p3 = load('T5_27_01_19.asvp');
p4 = load('T5_30_01_19.asvp');

Average_OBS_Depth = (1455+1496+1345+1336+1460+1239+1087+1370+1092+1025+1367+1157+916+1433+1475)/15;
%lower limit is at 15 m (source depth) and upper limit is at 1283m(average
%obs depth)

h = Average_OBS_Depth-15;
sum1 = 0;
sum2 = 0;
sum3 = 0;

for i = 22:1983 % for station XVT2
    
    g = (p1(i,2)+p1(i+1,2))/2;

    avgtdiff = (p1(i+1,1)-p1(i,1))/(g);

    sum1 = sum1+avgtdiff;
    clear avgtdiff g

end

v1 = h/sum1

clear i g avgtdiff summ

for i = 22:1817 % for station XVT5
    
    g = (p3(i,2)+p3(i+1,2))/2;

    avgtdiff = (p3(i+1,1)-p3(i,1))/(g);

    sum2 = sum2 + avgtdiff;
    clear avgtdiff g

end
v2 = h/sum2;

clear i g avgtdiff

for i = 22:1984 % for station XVT6
    
    g = (p4(i,2)+p4(i+1,2))/2;

    avgtdiff = (p4(i+1,1)-p4(i,1))/(g);

    sum3 = sum3+avgtdiff;
    clear avgtdiff
    
end

v3 = h/sum3;

v = (v1+v2+v3)/3;

figure(2)
plot(p1(22:1983,2),-1*p1(22:1983,1))
hold on
plot(p3(22:1818,2),-1*p3(22:1818,1))
plot(p4(22:1984,2),-1*p4(22:1984,1))
legend('XVT2','XVT5','XVT6')

clear i g avgtdiff summ Average_OBS_Depth h p1 p3 p4