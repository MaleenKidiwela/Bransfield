function [time, phase, status] = load_tlArrivalPickMatrix(fileStruc, traceMetaData, phase)
% Loads tlArrival times into a column matrix
%
% Usage 
%   [time, phase, status] = load_tlArrivalPickMatrix(fileStruc, traceMetaData, phase)
% Inputs
%   fileStruc   - tlArrival structure of filename of .mat with one
%   traceMetaData - Trace metadata with the following fields required
%                      station
%                      eventid
%                      ntrace*
%                      stationList*
%   phase         - Cell stucture of phase names (or string for one)
%                   Default option is all phases in tlArrival
% Outputs
%   time          - column matrix of times with one column per phase and
%                   entry per column per trace (NaN for no time)
%   phase         - Cell structure of phases names for each column
%   status        - Status of program
%                    1 - No file
%                    2 - Some traces have >1 time for at least one phase
%                        (none will be returned for such traces)


%% Process Inputs
% Get the tlArrival structure
if isstruct(fileStruc)
  eval(['tlArrival = ' fileStruc ';'])
else
  j = findstr(fileStruc,'.mat');
  if isempty(j)
    fileStruc = [fileStruc '.mat'];
  end
  if exist(fileStruc)==2
    tlArrival = load_tlArrival(fileStruc);
  else
    status = 1;
    time = [];
    return
  end
end

% Allow for old field 'type'
if ~isfield(tlArrival,'eventtype')
  tlArrival.eventtype = tlArrival.type;
  rmfield(tlArrival, 'type');
end

% Create missing arguments in traceMetaData
if ~isfield(traceMetaData,'ntrace')
  traceMetaData.ntrace = length(traceMetaData.eventid);
end
if ~isfield(traceMetaData,'stationList')
  traceMetaData.stationList = unique(traceMetaData.station);
end

% If phase is not specified get all phases in tlArrival
if nargin<3
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
    j = find(strcmp(traceMetaData.station, traceMetaData.stationList{is}));
    tlArrivalSubset = subset_tlArrival_tlPicker(tlArrival, phase{ip}, traceMetaData.stationList{is}, [], []);
    index = vector_indexmatch(tlArrivalSubset.eventid, traceMetaData.eventid(j));
    % >1 times for given phase
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