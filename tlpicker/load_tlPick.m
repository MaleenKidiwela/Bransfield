function [tlPick, status] = ... 
  load_tlPick(dirName, tMD, onePhase, channelSpecific, fieldList)
% Load picks for given traces and one phase from ASCII tlPick files into a tlPick structure
% 
% Usage 
%   [tlPick, status] = ... 
%           load_tlPick(dirName, tMD, onePhase, channelSpecific, fieldList);
% 
% Inputs
%   dirName         - Directory containing tlPick files
%   tMD             - Trace meta data for required traces
%                     The following fieldList are required in tMD
%                       station
%                       eventid
%                       channel - only if channelSpecific == True
%                     The following fieldList are created if not present
%                       ntrace
%                       stationList
%   onePhase      -   String or cell with a single phase name
%   channelSpecific - Pick channel needs to match tMD.channel
%   fieldList       - Fields to return in tlPick in a cell array
%                     Options are
%                        all
%                        station
%                        channel
%                        eventid
%                        phase
%                        time
%                        unc
%                        filtLim0
%                        filtLim1
%                        filtOrder
%                        filtZeroPhase
%                        scale
%                        user
%                        lddate
%                        use
%                        comment
%                        stationList
%                        phaseList
%                        picked
%                        updated
%
% Outputs
%   tlPick - Structure with fields requested by fieldNames
%   status - Status of execution
%            0 - okay
%            1 - tlPick directory does not exist
% 4/22/11 - Force station and stationList fields of tMD to be a cell structure

status = 0;

%% Process inputs
if ~iscell(tMD.station)
  tMD.station = {tMD.station};
end
if iscell(onePhase)
  onePhase = cell2mat(onePhase);
end
if nargin<4
  channelSpecific = false;
end
if nargin<5
  fieldList = {'all'};
end
if exist(dirName)~=7 && ~isempty(dirName)
  warning(['save_tlPick: Directory ' dirName ' does not exist']);
  status = 1;
  return
elseif ~isempty(dirName)
  if dirName(end)~='/'
    dirName = [dirName '/'];
  end
end
if ~isfield(tMD,'ntrace')
  tMD.ntrace = length(tMD.eventid);
end
if ~isfield(tMD,'channel') && ~channelSpecific
  tMD.channel = zeros(tMD.ntrace,1);
end
if ~isfield(tMD,'stationList')
  tMD.stationList = unique(tMD.station);
elseif ~iscell(tMD.stationList)
  tMD.stationList = {tMD.stationList};
end
if channelSpecific
  channelList = unique(tMD.channel);
else
  channelList = 0;
end

%% Create pick fieldList
if any(strcmpi(fieldList,'station')) || any(strcmpi(fieldList,'all'));
  tlPick.station = tMD.station;
end
if any(strcmpi(fieldList,'npick')) || any(strcmpi(fieldList,'all'));
  tlPick.npick = tMD.ntrace;
end
if any(strcmpi(fieldList,'channel')) || any(strcmpi(fieldList,'all'));
  tlPick.channel = tMD.channel;
end
if any(strcmpi(fieldList,'eventid')) || any(strcmpi(fieldList,'all'));
  tlPick.eventid = tMD.eventid;
end
if any(strcmpi(fieldList,'phase')) || any(strcmpi(fieldList,'all'));
  tlPick.phase = cellstr(repmat(onePhase,tMD.ntrace,1));
end
if any(strcmpi(fieldList,'time')) || any(strcmpi(fieldList,'all'));
  tlPick.time = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'unc')) || any(strcmpi(fieldList,'all'));
  tlPick.unc = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'filtLim0')) || any(strcmpi(fieldList,'all'));
  tlPick.filtLim0 = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'filtLim1')) || any(strcmpi(fieldList,'all'));
  tlPick.filtLim1 = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'filtOrder')) || any(strcmpi(fieldList,'all'));
  tlPick.filtOrder = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'filtZeroPhase')) || any(strcmpi(fieldList,'all'));
  tlPick.filtZeroPhase = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'scale')) || any(strcmpi(fieldList,'all'));
  tlPick.scale = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'user')) || any(strcmpi(fieldList,'all'));
  tlPick.user = cellstr(repmat(' ',tMD.ntrace,1));
end
if any(strcmpi(fieldList,'lddate')) || any(strcmpi(fieldList,'all'));
  tlPick.lddate = NaN(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'use')) || any(strcmpi(fieldList,'all'));
  tlPick.use = true(tMD.ntrace,1);
end
if any(strcmpi(fieldList,'comment')) || any(strcmpi(fieldList,'all'));
  tlPick.comment = cellstr(repmat(' ',tMD.ntrace,1));
end

%% Loop through getting tlPick
for currentStation = tMD.stationList(:)'

  [station, channel, eventid, phase, time, unc, use,  ...
  filtLim0, filtLim1, filtOrder, filtZeroPhase, scale, user, ...
  lddate, comment, statusRead] = read_tlPick(dirName, currentStation, onePhase, {'all'});

  if ~statusRead

    for currentChannel = channelList

      if channelSpecific
        i2 = find(strcmp(tMD.station,currentStation) & tMD.channel==currentChannel);
        i1 = find(strcmp(station,currentStation) & channel==currentChannel);
      else
        i2 = find(strcmp(tMD.station,currentStation));
        i1 = find(strcmp(station,currentStation));
      end
      j = vector_indexmatch(eventid(i1),tMD.eventid(i2));
      if any(j<0);
        warning(['load_tlPick: Found duplicate tlPick']);
        error('Aborting')
      end
      i2 = i2(j>0);
      j = j(j>0);

      if any(strcmpi(fieldList,'channel')) || any(strcmpi(fieldList,'all'));
        tlPick.channel(i2) = channel(i1(j));
      end
      if any(strcmpi(fieldList,'phase')) || any(strcmpi(fieldList,'all'));
        tlPick.phase(i2) = phase(i1(j));
      end
      if any(strcmpi(fieldList,'time')) || any(strcmpi(fieldList,'all'));
        tlPick.time(i2) = time(i1(j));
      end
      if any(strcmpi(fieldList,'unc')) || any(strcmpi(fieldList,'all'));
        tlPick.unc(i2) = unc(i1(j));
      end
      if any(strcmpi(fieldList,'filtLim0')) || any(strcmpi(fieldList,'all'));
        tlPick.filtLim0(i2) = filtLim0(i1(j));
      end
      if any(strcmpi(fieldList,'filtLim1')) || any(strcmpi(fieldList,'all'));
        tlPick.filtLim1(i2) = filtLim1(i1(j));
      end
      if any(strcmpi(fieldList,'filtOrder')) || any(strcmpi(fieldList,'all'));
        tlPick.filtOrder(i2) = filtOrder(i1(j));
      end
      if any(strcmpi(fieldList,'filtZeroPhase')) || any(strcmpi(fieldList,'all'));
        tlPick.filtZeroPhase(i2) = filtZeroPhase(i1(j));
      end
      if any(strcmpi(fieldList,'scale')) || any(strcmpi(fieldList,'all'));
        tlPick.scale(i2) = scale(i1(j));
      end
      if any(strcmpi(fieldList,'user')) || any(strcmpi(fieldList,'all'));
        tlPick.user(i2) = user(i1(j));
      end
      if any(strcmpi(fieldList,'lddate')) || any(strcmpi(fieldList,'all'));
        tlPick.lddate(i2) = lddate(i1(j));
      end
      if any(strcmpi(fieldList,'use')) || any(strcmpi(fieldList,'all'));
        tlPick.use(i2) = use(i1(j));
      end
      if any(strcmpi(fieldList,'comment')) || any(strcmpi(fieldList,'all'));
        tlPick.comment(i2) = comment(i1(j));
      end

    end %currentChannel

  end

end %currentStation

%% Final assignments
if any(strcmpi(fieldList,'stationList')) || any(strcmpi(fieldList,'all'));
  tlPick.stationList = tMD.stationList;
end
if any(strcmpi(fieldList,'phaseList')) || any(strcmpi(fieldList,'all'));
  tlPick.phaseList = {onePhase};
end
if any(strcmpi(fieldList,'picked')) || any(strcmpi(fieldList,'all'));
  tlPick.picked = ~isnan(tlPick.time);
end
if any(strcmpi(fieldList,'updated')) || any(strcmpi(fieldList,'all'));
  tlPick.updated = false(tlPick.npick,1);
end
