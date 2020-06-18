function srStationRotation = load_srStationRotation(theFile,srStation);

if ~exist(theFile)
  warning(['Rotation file ' theFile ' does not exist'])
  for i = 1:srStation.nsta
    srStationRotation(i).name = srStation.name(i);
    srStationRotation(i).use = false;
    srStationRotation(i).rot = [];
    srStationRotation(i).scale = [];
    % Optional
    srStationRotation(i).processed = false;
    srStationRotation(i).meanMisfit = [];
    srStationRotation(i).medianMisfit = [];
    srStationRotation(i).rotAlt = [];
    srStationRotation(i).scaleAlt = [];
    srStationRotation(i).meanMisfitAlt = [];
    srStationRotation(i).medianMisfitAlt = [];
    srStationRotation(i).segy = [];
    srStationRotation(i).filter = [];
    srStationRotation(i).static = [];
    srStationRotation(i).active = [];
    srStationRotation(i).orientHoriz = [];
    srStationRotation(i).user = [];
  end  
else
  eval(['load ' theFile])
  index = zeros(srStation.nsta,1);
  for i = 1:length(srStationRotation)
    j = find(strcmpi(srStation.name,srStationRotation(i).name));
    if length(j)==1
      index(i) = j;
    end
  end
  if ~all(index)
    error('srStationRotation does not match srStation')
  else
    srStationRotation = srStationRotation(index);
  end
end