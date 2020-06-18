function upicker2tlPick(user);
%Function to convert upicker files to tlpick format
% Usage
%   upicker2tlPick(user)
% 
% Inputs
%   user - Username which defaults to 'unknown'
%
%  upicker filenames have format *_obs##*_@@.mat where ## is instrument number
%     and @@ is phase name
%  All files in current directory are converted
%
% Added code from EEH to deal with "_" in the phase name (3/31/11)
% Converted to a function so user can be specified (4/15/11)
% Converts single digits stations to single digit names (e.g., '05' to '5')
% (6/15/11)

if nargin<1
  user = 'unknown';
end
dirList = dir;
tic

for idir = 1:length(dirList)
  name = dirList(idir).name;
  i1 = strfind(name,'obs');
  i2 = strfind(name,'.mat');
  if length(i1)==1 && length(i2)==1
    disp(['upicker2tlPick: Converting file ' name]);
    name = name(i1+3:end);
    i2 = i2-i1-2;
    i3 = findstr(name,'_');
    station = name(1:i3(1)-1);
    if station(1) == '0';
      station = station(2:end);
    end
    if length(i3) == 2,
        phase = name(i3(end)+1:i2-1);
    else   % Case where the phase name has an underscore in it EEH
        phase = [name(i3(end-1)+1:i3(end)-1) name(i3(end)+1:i2-1)];
        disp(['Changed code to change phase name to ' phase])
    end    
    clear pick
    eval(['load ' dirList(idir).name]);
    npick = length(tpick);
    pick.station = repmat({station},npick,1);
    pick.channel = zeros(npick,1);
    pick.eventid = idpick(:);
    pick.phase = repmat({phase},npick,1);
    pick.time = tpick(:);
    pick.unc = errpick(:);
    pick.user = repmat({user},npick,1);
    status = save_tlPick('.',pick);
  end
end
toc
    
    
  
