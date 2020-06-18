function tlOutput2tlPick(srEvent, srStation, tlOutput, iteration, phaseIn, phaseOut, tlPickDir)
% Saves predicted times from a tomolab tlOutput directory to tlPick structures 
%
% Usage
%   tlOutput2tlPick(srEvent, srStation, tlOutput, iteration, phaseIn, phaseOut, tlPickDir)
%
% Inputs
%   srEvent   - srEvent structure with derived x-y variable (use load_srEvent)
%   srStation - srStation structure with derived x-y variables (use load_srStation)
%   tlOuput   - directory with tomolab output
%   iteration - iteration of a tlMisfit (1 = forward model)
%   phaseIn   - cell array of phases in tlPick files
%               (empty implies all phases in tlMisfit structure)
%   phaseOut  - Corresponding phase names in tlArrival structure 
%               (empty implies phaseOut = phaseIn)
%   tlPickDir - Directory in which to save tlPick files (Default = '.')
%

global debugPicking
debugPicking = false;

if nargin<6
  phaseOut = [];
end
if nargin<7
  tlPickDir = '.';
end
if isempty(tlPickDir)
  tlPickDir = '.';
end
if tlPickDir(end) ~= '/'
  tlPickDir = [tlPickDir '/'];
end

% Create traceMetaData input for load_tlOutputPickMatrix
tmd.station = repmat(srStation.name', srEvent.nevt, 1);
tmd.station = tmd.station(:);
tmd.eventid = repmat(srEvent.id, srStation.nsta, 1);

% Get times and phases
[time, phaseIn, status] = load_tlOutputPickMatrix(tlOutput, tmd, iteration, phaseIn);

if isempty(phaseOut) 
  phaseOut = phaseIn;
end

comment = [tlOutput ' Iter=' int2str(iteration)];
comment = comment(1:min(80,length(comment)));
for iPhase = 1:length(phaseOut)
  for iStation = 1:srStation.nsta
    i = strcmp(tmd.station,srStation.name{iStation});
    if ~all(isnan(time(i,iPhase)));
      tlPick.time = time(i,iPhase);
      tlPick.station = tmd.station(i);
      tlPick.eventid = tmd.eventid(i);
      tlPick.channel = zeros(length(tlPick.time),1);
      tlPick.phase = repmat(phaseOut(iPhase),length(tlPick.time),1);
      tlPick.unc = zeros(length(tlPick.time),1);
      tlPick.comment = repmat({comment},length(tlPick.time),1);
      
      filename = [tlPickDir 'tlPick_' srStation.name{iStation} '_' phaseOut{iPhase} '.dat'];
      if exist(filename)==2
        eval(['delete ' filename]);
      end
      status = save_tlPick(tlPickDir, tlPick, false, true);
      disp(['Writing tlPick file ' filename])
    end 
  end
end
