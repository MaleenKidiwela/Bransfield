function srElevation = convert_grd2srElevation(grdFilename,fudge)
% Converts a GRD bathymetry file to Stringray elevation format 
%
% Usage
%   srElevation = convert_grd2srElevation(grdFilename)
% Inputs
%   grdFilename - GRD file name
% Outputs
%   srElevation - Stingray elevation structure
%   fudge - logical to remove N-S lines of NaN's in Endeavour_prelim.grd 
%
% Uses GMT commands
grdFilename = 'GMRTv3_7_20200519topo.grd';

status = unix(['/Users/earthnote/Desktop/Brainsfield/' 'GMRTv3_7_20200519topo.grd' '> temp.txt'])

if status
  error(['convert_grd2srElevation(' grdFilename ') could not execute grd2xyz'])
end

load temp.txt -ascii
!rm temp.txt
x_min = min(temp(:,1));
x_max = max(temp(:,1));
nx = find(diff(temp(:,1))<0,1,'first');
x_inc = (x_max - x_min) / (nx-1);
y_min = min(temp(:,2));
y_max = max(temp(:,2));
ny = size(temp,1) / nx;
y_inc = (y_max - y_min) / (ny-1);

data = reshape(temp(:,3),nx,ny);
% Make sure latitude is increasing with second index
if temp(1,2)-temp(ny+1,2)>0
  data = fliplr(data);
end

if fudge
  for i = [1484 1493 1503 1513 1523]
    data(i,:) = (data(i-1,:) + data(i+1,:))/2;
  end
end

srElevation.header = [x_min x_max y_min y_max x_inc y_inc nx ny];
srElevation.data = data;