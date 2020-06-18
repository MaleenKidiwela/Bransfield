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
  filtLim0 = 0;
  nfiltLim0 = 1;
else
  filtLim0 = pick.filtLim0;
  nfiltLim0 = npick;
end
if ~isfield(pick,'filtLim1')
  filtLim1 = zeros(npick,1);
else
  filtLim1 = pick.filtLim1;
end
if ~isfield(pick,'filtOrder')
  filtOrder = zeros(npick,1);
else
  filtOrder = pick.filtOrder;
end
if ~isfield(pick,'filtZeroPhase')
  filtZeroPhase = zeros(npick,1);
else
  filtZeroPhase = pick.filtZeroPhase;
end
if ~isfield(pick,'scale')
  scale = zeros(npick,1);
else
  scale = pick.scale;
end
if ~isfield(pick,'user')
  user = repmat({'unknown'},npick,1);
else
  user = pick.user;
end
if ~isfield(pick,'lddate')
  c = clock;
  lddate = repmat(date2secnds(c(1),julday(c(2),c(3),c(1)),c(4),c(5),c(6)),npick,1);
else
  lddate = pick.lddate;
end
if ~isfield(pick,'comment')
  comment = repmat({'None'},npick,1);
else
  comment = pick.comment;
end

if ~isfield(pick,'stationList')
  stationList = unique(pick.station);
else
  stationList = pick.stationList;
end
if ~isfield(pick,'phaseList')
  phaseList = unique(pick.phase);
else
  phaseList = pick.phaseList;  
end
toc

for station = stationList
  for phase = phaseList
    iout = strcmp(pick.phase,phase) & strcmp(pick.phase,phase);
    if any(iout)
      iout = find(iout);
      fileName = [dirName 'tlPick_' station{1} '_' phase{1} '.dat'];
      if ~exist(fileName)
        fid = fopen(fileName,'w');
        append = 0;
      else
        % APPENDING HERE  
      end
      toc
      
      c = zeros(1,212*length(iout),'uint8');
      j = 0;
      for i = iout(:)'
%         fprintf(fid,'%6s%2i%8i%8s%10.4f%7.4f%5.1f%5.1f%2i%2i%10.2e%10s%17.5f%40s%80s', ...
%                 pick.station{i},pick.channel(i),pick.eventid(i),pick.phase{i}, ...
%                 pick.time(i),pick.unc(i),filtLim0(i),filtLim0(i), ...
%                 filtOrder(i),filtZeroPhase(i),scale(i),user{i}, ...
%                 lddate(i), ' ', comment{i});
        c(j+1:j+212) = sprintf('%6s%2i%8i%8s%10.4f%7.4f%5.1f%5.1f%2i%2i%10.2e%10s%17.5f%40s%80s', ...
                pick.station{i},pick.channel(i),pick.eventid(i),pick.phase{i}, ...
                pick.time(i),pick.unc(i),filtLim0(i),filtLim0(i), ...
                filtOrder(i),filtZeroPhase(i),scale(i),user{i}, ...
                lddate(i), ' ', comment{i});
        j = j + 212;
      end
      fwrite(fid,c,'uint8');
      fclose(fid);
    end
  end
end
 

  

  