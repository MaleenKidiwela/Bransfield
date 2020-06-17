% first run ---- gmt xyz2grd Orca15m.txt -R-62.6341/-62.2814/-58.7208/-57.9002 I-0.000549 -Gbravoseis
% then run ----- gmt grd2xyz bravoseis.grd > bravoseis.xyz 
% on gmt before executing this


load bravoseis.xyz -ascii
% !rm x.xyz
x_min = min(bravoseis(:,1));
x_max = max(bravoseis(:,1));
nx = find(diff(bravoseis(:,1))<0,1,'first');
x_inc = (x_max - x_min) / (nx-1);
y_min = min(bravoseis(:,2));
y_max = max(bravoseis(:,2));
ny = size(bravoseis,1) / nx;
y_inc = (y_max - y_min) / (ny-1);


data = reshape(bravoseis(:,3),nx,ny);

% Make sure latitude is increasing with second index
if bravoseis(1,2)-bravoseis(ny-1,2)>0 % ny+1 originally here brings an error 
  data = fliplr(data);
end


srElevation.header = [x_min x_max y_min y_max x_inc y_inc nx ny];
srElevation.data = data;