function [station, channel, eventid, phase, time, unc, use,  ...
  filtLim0, filtLim1, filtOrder, filtZeroPhase, scale, user, ...
  lddate, comment, status] = read_tlPick(dirName, stationInput, onePhase, fieldList)
% Read one ASCII tlPick_STATION_PHASE.dat file
%
% Usage
%  [station, channel, eventid, phase, time, unc,  ...
%    filtLim0, filtLim1, filtOrder, filtZeroPhase, scale, user, ...
%    lddate, comment, status] = read_tlPick(dirName, stationInput, onePhase, fieldList)
%
% Inputs
%   dirName      - Directory containing tlPick files
%   stationInput - Station string for a single station
%   onePhase     - Phase string for a single phase
%   fieldList    - Cell array with fields to read ( Use {'all'} for all fields)
%                  Default is {'all'}
%                  Fields not specified are returned as an empty variable
%
% Outputs 
%    Columns in pick file (empty if not requested in fieldList)
%      station
%      channel
%      eventid
%      phase
%      time
%      unc
%      use
%      filtLim0
%      filtLim1
%      filtOrder
%      filtZeroPhase
%      scale
%      user
%      lddate
%      comment
%    status - Status of execution
%             0 - Okay
%             1 - tlPick directory does not exist
%             2 - tlPick file does not exist
%             3 - tlPick file exists but is empty
%
% Outputs are not put into a scalar structure of vectors and cell vectors
% since this really slows things down

% Parameters
lengthRecord = 213;

% Process inputs
if ~isempty(dirName)
  if dirName(end) == '/'
    dirName = dirName(1:end-1);
  end
else
  dirName = '.';
end
if exist(dirName)~=7
  warning(['read_tlPick: Directory ' dirName ' does not exist']);
  status = 1;
  return
end

if iscell(stationInput)
  stationInput = cell2mat(stationInput);
end
if iscell(onePhase)
  onePhase = cell2mat(onePhase);
end
fileName = [dirName '/tlPick_' stationInput '_' onePhase '.dat'];
if nargin<4; fieldList = {'all'}; end

% Default outputs
status = 0;
station = [];
channel = [];
eventid = [];
phase = [];
time = [];
unc = []; 
use = [];
filtLim0 = [];
filtLim1 = [];
filtOrder = [];
filtZeroPhase = [];
scale = [];
user = [];
lddate = [];
comment = [];

% Read data 
if exist(fileName)~=2
  warning(['read_tlPick - tlPick file ' fileName ' does not exist'])
  status = 2;
  return
else
  fid = fopen(fileName,'r');
  [c,n] = fread(fid,inf,'uint8');
  m = n/lengthRecord;
  c = reshape(c,lengthRecord,m);
%   string = char(c(1:1,:));
%   notDeleted = ~sscanf(string(:)','%1i');
%   byte0 = (find(notDeleted)-1)*lengthRecords;
  if isempty(c)
    status = 3;
  else
    if any(strcmpi(fieldList,'station')) || any(strcmpi(fieldList,'all'))
      station = deblank_fb(cellstr(char(c(1:6,:))'));
    end
    if any(strcmpi(fieldList,'channel')) || any(strcmpi(fieldList,'all'))
      string = char(c(7:8,:));
      channel = sscanf(string(:),'%2i');  
    end
    if any(strcmpi(fieldList,'eventid')) || any(strcmpi(fieldList,'all'))
      string = char(c(9:16,:));
      eventid = sscanf(string,'%8i');  
    end
    if any(strcmpi(fieldList,'phase')) || any(strcmpi(fieldList,'all'))
      phase = deblank_fb(cellstr(char(c(17:24,:))'));
    end
    if any(strcmpi(fieldList,'time')) || any(strcmpi(fieldList,'all'))
      string = char(c(25:34,:));
      time = sscanf(string(:)','%10f');  
    end
    if any(strcmpi(fieldList,'unc')) || any(strcmpi(fieldList,'all'))
      string = char(c(35:41,:));
      unc = sscanf(string,'%7f');  
    end
    if any(strcmpi(fieldList,'filtLim0')) || any(strcmpi(fieldList,'all'))
      string = char(c(42:46,:));
      filtLim0 = sscanf(string,'%5f');  
    end
    if any(strcmpi(fieldList,'filtLim1')) || any(strcmpi(fieldList,'all'))
      string = char(c(47:51,:));
      filtLim1 = sscanf(string,'%5f');  
    end
    if any(strcmpi(fieldList,'filtOrder')) || any(strcmpi(fieldList,'all'))
      string = char(c(52:53,:));
      filtOrder = sscanf(string,'%2i');  
    end
    if any(strcmpi(fieldList,'filtZeroPhase')) || any(strcmpi(fieldList,'all'))
      string = char(c(54:55,:));
      filtZeroPhase = sscanf(string,'%2i');  
    end
    if any(strcmpi(fieldList,'scale')) || any(strcmpi(fieldList,'all'))
      string = char(c(56:65,:));
      scale = sscanf(string,'%10e');  
    end
    if any(strcmpi(fieldList,'user')) || any(strcmpi(fieldList,'all'))
      user = deblank_fb(cellstr(char(c(66:75,:))'));
    end
    if any(strcmpi(fieldList,'lddate')) || any(strcmpi(fieldList,'all'))
      string = char(c(76:92,:));
      lddate = sscanf(string,'%17f');  
    end
    if any(strcmpi(fieldList,'use')) || any(strcmpi(fieldList,'all'))
      string = char(c(94:94,:));
      string(abs(string)==32) = '1';            % Default is use for backward compatability
      use = logical(sscanf(string,'%1i'));  
    end
    if any(strcmpi(fieldList,'comment')) || any(strcmpi(fieldList,'all'))
      comment = deblank(cellstr(char(c(133:212,:))'));
    end
  end
  fclose(fid);
end
  
