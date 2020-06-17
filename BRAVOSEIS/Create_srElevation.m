% run gmt grd2xyz GMRTv3_7_20200519topo.grd > x.xyz on gmt before executing
% this

load x.xyz -ascii
% !rm x.xyz
x_min = min(x(:,1));
x_max = max(x(:,1));
nx = find(diff(x(:,1))<0,1,'first');
x_inc = (x_max - x_min) / (nx-1);
y_min = min(x(:,2));
y_max = max(x(:,2));
ny = size(x,1) / nx;
y_inc = (y_max - y_min) / (ny-1);

%%
data = reshape(x(:,3),nx,ny);

% Make sure latitude is increasing with second index
if x(1,2)-x(ny-1,2)>0 % ny+1 originally here brings an error 
  data = fliplr(data);
end


srElevation.header = [x_min x_max y_min y_max x_inc y_inc nx ny];
srElevation.data = data;