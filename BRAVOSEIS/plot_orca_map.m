% load Orca15m.txt -ascii % make sure to unzip the .zip file
% x = [Orca15m(:,1)];
% y = [Orca15m(:,2)];
% z = [Orca15m(:,3)];
% figure(1)
% stem3(x, y, z)
% grid on
% xv = linspace(min(x), max(x), 100);
% yv = linspace(min(y), max(y), 100);
% [X,Y] = meshgrid(xv, yv);
% Z = griddata(x,y,z,X,Y);
% figure(2)
% surf(X, Y, Z);
% grid on
% set(gca, 'ZLim',[0 1800])
% shading interp

%plot below
create_bravoseis_grd_xyz
a= isnan(bravoseis(:,3));
tab = readtable('orca_tomo_shotfile_final.txt');
plot(bravoseis(a,1),bravoseis(a,2),'.')
hold on
plot(tab.sourceLat,tab.sourceLon,'.')

