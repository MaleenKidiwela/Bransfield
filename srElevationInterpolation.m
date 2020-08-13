%interpolating NaNs within the srElevation file made from 15m orca txt
%bathymetric data

load('srElevation_orca.mat')
load('srStation_orca.mat')
load('srEvent_orca.mat')

[X,Y,z_grd] = grdread2('GMRTv3_7_20200708topo.grd');
%[yy,xx,zz]= grdread2('bravoseis.grd');
%zz=zz';
Y = Y'
%yy = yy'
figure
surf(srElevation.LON,srElevation.LAT,srElevation.data,'edgecolor', 'none')
hold on
surf(X,Y,z_grd,'edgecolor', 'none')


for i = 1:length(Y)
    
    x_grd(i,:) = X;

end

for j = 1:length(X)
    
    y_grd(:,j) = Y;

end


% 
% for i = 1:length(yy)
%     
%     srElevation.LON(i,:) = xx;
% 
% end
% 
% for j = 1:length(xx)
%     
%     srElevation.LAT(:,j) = yy;
% 
% end
% 
% srElevation.data = zz;

lenElev = size(srElevation.data);
intp = interp2(x_grd,y_grd,z_grd,srElevation.LON,srElevation.LAT);
for x = 1: lenElev(1)
    
    for y = 1: lenElev(2)
        
        if isnan(srElevation.data(x,y))
            
            
            srElevation.data(x,y)= intp(x,y);
            
        end
        
    
    end
    (x/lenElev(1))*100
end

figure
surf(srElevation.LON,srElevation.LAT,srElevation.data,'edgecolor', 'none')
hold on
scatter3(srStation.longitude,srStation.latitude,srStation.elevation*1000,'MarkerFaceColor','red')
plot(srEvent.longitude,srEvent.latitude,'r.')

