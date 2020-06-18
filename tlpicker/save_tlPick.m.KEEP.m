function status = save_tlPick(dirName, pick, channelSpecific)

status = 0;

if nargin<3
  channelSpecific = 0;
end

if exist(dirName)~=7
  disp(['save_tlPick: Directory ' dirName ' does not exist']);
  status = 1;
  return
else
  if dirName(end)~='/'
    dirName = [dirName '/'];
  end
end
if ~all(isfield(pick,{'station','eventid','phase','time','unc','channel'}))
  disp(['save_tlPick: At least one required pick fields does not exist']);
  status = 2;
  return
end

npick = length(pick.time);
if ~isfield(pick,'filtLim0')
  pick.filtLim0 = zeros(npick,1);
end
if ~isfield(pick,'filtLim1')
  pick.filtLim1 = zeros(npick,1);
end
if ~isfield(pick,'filtOrder')
  pick.filtOrder = zeros(npick,1);
end
if ~isfield(pick,'filtZeroPhase')
  pick.filtZeroPhase = zeros(npick,1);
end
if ~isfield(pick,'scale')
  pick.scale = zeros(npick,1);
end
if ~isfield(pick,'user')
  pick.user = repmat({'unknown'},npick,1);
end
if ~isfield(pick,'lddate')
  c = clock;
  pick.lddate = repmat(date2secnds(c(1),julday(c(2),c(3),c(1)),c(4),c(5),c(6)),npick,1);
end
if ~isfield(pick,'comment')
  pick.comment = repmat({'None'},npick,1);
end

if ~isfield(pick,'stationList')
  pick.stationList = unique(pick.station);
end
if ~isfield(pick,'phaseList')
  pick.phaseList = unique(pick.phase);
end

for station = pick.stationList
  for phase = pick.phaseList
    iout = strcmp(pick.phase,phase) & strcmp(pick.phase,phase);
    if any(iout)
      iout = find(iout);
      fileName = [dirName 'tlPick_' station{1} '_' phase{1} '.dat'];
      if ~exist(fileName)
        fid = fopen(fileName,'wt');
        append = 0;
      else
        % APPENDING HERE  
      end
      
      for i = iout(:)'
        fprintf(fid,'%6s%2i%8i%8s%10.4f%7.4f%5.1f%5.1f%2i%2i%10.2e%10s%17.5f%40s%80s', ...
                pick.station{i},pick.channel(i),pick.eventid(i),pick.phase{i}, ...
                pick.time(i),pick.unc(i),pick.filtLim0(i),pick.filtLim0(i), ...
                pick.filtOrder(1),pick.filtZeroPhase(i),pick.scale(i),pick.user{i}, ...
                pick.lddate, ' ', pick.comment{i});
      end
      fclose(fid);
    end
  end
end
 

  

  