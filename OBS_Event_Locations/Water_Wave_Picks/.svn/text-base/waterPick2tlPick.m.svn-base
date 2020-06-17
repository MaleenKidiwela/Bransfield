% Script to convert Anne Wells' water wave picks to tlPick files
% Need to set input and output directories
%
% Will need to deal with timing corrections for a couple of OBSs

waterPickDir = '/Users/matlab/OBS_Event_Locations/trunk/Water_Wave_Picks';
tlPickDir = '/users/matlab/etomo_picks/Picks';

if waterPickDir(end)~='/' 
  waterPickDir = [waterPickDir '/'];
end

for i=1:68
  % Convert Annes integer station names to string names
  if i<=64
    station = int2str(i);
  elseif i==65
    station = '15A';
  elseif i==66
    station = '15B';
  elseif i==67
    station = '15C';    
  elseif i==68
    station = '46A';
  end
  
  % Load her picks and convert
  eval(['load ' waterPickDir 'waterPick_OBS_' int2str(i) '.mat']);
  w = waterPick(i);
  keep = ~isnan(w.artime);
  clear tlPick
  tlPick.npick = sum(keep);
  tlPick.station = repmat({station},tlPick.npick,1);
  tlPick.channel = ones(tlPick.npick,1);
  tlPick.phase = repmat({'W'},tlPick.npick,1);
  tlPick.time = w.artime(keep);
  tlPick.eventid = w.shot(keep);
  tlPick.unc =  w.sdv(keep);
  status = save_tlPick(tlPickDir, tlPick, false, true);
  
end 
  
