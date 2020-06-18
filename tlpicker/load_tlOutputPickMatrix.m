function [time, phase, status] = load_tlOutputPickMatrix(tlOutput, traceMetaData, iteration, phase)
% Loads tlMisfit predictions from a tlOutput directory into a column matrix
%
% Usage 
%   [time, phase, status] = load_tlOutputPickMatrix(tlOutput,traceMetaData, iteration, phase)
% Inputs
%   tlOutput      - Directory with tomolab outputs
%   traceMetaData - Trace metadata with the following fields required
%                      station
%                      eventid
%                   and the following that will be derived if not present
%                      ntrace
%                      stationList
%   iteration     - Inversion iteration number for the tlMisfit file
%                   (default 1)
%   phase         - Cell stucture of phase names (or string for one)
%                   Default option is all phases in tlArrival
% Outputs
%   time          - column matrix of times with one column per phase and
%                   entry per column per trace (NaN for no time)
%   phase         - Cell structure of phases names for each column (default
%                   is all the phases in the file)
%   status        - Status of program
%                    1 - No file
%                    2 - Some traces have >1 time for at least one phase
%                        (none will be returned for such traces)

%% Process Inputs
if exist(tlOutput)~=7
  warning(['load_tlOutputPickMatrix: tlOutput directory ' tlOutput ' does not exist'])
  time = [];
  phase = [];
  status = 1;
  return
end
if nargin<3
  interation = 1;
end
if isempty(iteration)
  iteration = 1;
end
if nargin<4
  phase = [];
end
% Create missing arguments in traceMetaData
if ~isfield(traceMetaData,'ntrace')
  traceMetaData.ntrace = length(traceMetaData.eventid);
end
if ~isfield(traceMetaData,'stationList')
  traceMetaData.stationList = unique(traceMetaData.station);
end


%% Load necessary structures
fileName = [tlOutput '/tlControl.mat'];
load(fileName)
fileName = [tlOutput '/tlMisfit_it' int2str(iteration) '.mat'];
load(fileName)
fileName = tlControl.files.Arrival;
tlArrival = load_tlArrival(fileName);

% Use pointer from tlMisfit to create corresponding tlArrival
% Allow for old field 'type'
tlArrival.time      = tlMisfit.ttime;
if isfield(tlArrival,'eventtype')
  tlArrival.eventtype = tlArrival.eventtype(tlMisfit.ptr);
else
  tlArrival.eventtype = tlArrival.type(tlMisfit.ptr);
  rmfield(tlArrival, 'type');
end
tlArrival.eventid   = tlArrival.eventid(tlMisfit.ptr);
tlArrival.phase     = tlArrival.phase(tlMisfit.ptr);
tlArrival.error     = tlArrival.error(tlMisfit.ptr);
tlArrival.station   = tlArrival.station(tlMisfit.ptr);

% If phase is not specified get all phases in tlArrival
if isempty(phase)
  phase = tlArrival.phaselist;
else
  if iscell(phase)
    if isempty(deblank(phase{1}))
      phase = tlArrival.phaselist;
    end
  else
    if isempty(deblank(phase))
      phase = tlArrival.phaselist;
    else
      phase = {phase};
    end
  end
end

%% Create outputs
% Preallocate matrix of times
time = NaN(traceMetaData.ntrace, length(phase));

% Loop through phases and stations getting times
status = 0;
for ip = 1:length(phase)
  for is = 1:length(traceMetaData.stationList);
    j = find(strcmp(traceMetaData.station, traceMetaData.stationList(is)));
    tlArrivalSubset = subset_tlArrival_tlPicker(tlArrival, phase{ip}, traceMetaData.stationList{is}, [], []);
    index = vector_indexmatch(tlArrivalSubset.eventid, traceMetaData.eventid(j));
    % >1 times for a given phase
    if any(index<0)
      status = 2;
    end
%     % 0 times for a given phase
%     if any(index==0) && ~rem(status,2)
%       status = status + 1;
%     end
    % Get times when there is only 1 
    k = find(index>0);
    time(j(k), ip) = tlArrivalSubset.time(index(k));
  end
end  
