function [pick, status] = get_pickMatrix(DFS, phaseList, chanIter, tMD)
% Create a column matrix of pick times for given traces and phases
%
% Works with picks in (1) tlPick ASCII files, 
%                     (2) tlArrival file/structure
%                  or (3) tlMisfit for an inversion directory
%
% Usage
%   [pick, status] = get_pickMatrix(DFS, phaseList, chanIter, tMD)
%
% Inputs
%   DFS       - Directory name for tlPick_STATION_PHASE.dat ASCII files
%                or Matlab file name / structure name for tlArrival
%                or directory name for inversion with tlMisfit
%   phaseList - String with one phase name or a cell array of Phase names 
%                corresponding to columns in pick
%                A blank entry gives no output for tlPick ASCII format but
%                all phases in a tlArrival or tlMisfit structure
%   chanIter  - For tlPick ASCII files this specifies how the pick channel
%               is treated with options as follows
%               0, NaN, Empty - Ignore the pick channel
%               <0            - Match pick channel to trace channel
%               >0            - Match picks for this channel to traces
%                               irrespective of trace channel
%               For the tlMisfit it is the inversion iteration number 
%                   of the tlMisfit file (default 1)
%   tMD       -  Trace metadata structure 
%                The following fields required
%                  station
%                  eventid
%                  channel for tlPick ASCII files and chanIter<0
%                The following fields will be created if they do not exist 
%                  ntrace
%                  stationList
%
% Outputs
%   pick   - Matrix of pick times with one column per phase and one row per
%            trace.  NaNs denote no pick 
%   status - Status of execution
%              0 - Okay
%              1 - DFS is an empty string
%              2 - DFS directory/file does not exist
%              3 - More than one pick for at least one trace 
%                 (tlArrival & tlMisfit only
%              4 - phase is empty (tlPick option only)

status = 0;
pick = [];

% Return if DFS is empty string
if isempty(DFS)
  status = 1;
  return
end

% Get inputs into required format
if iscell(DFS)
  DFS = cell2mat(DFS);
end
if ~iscell(phaseList)
  phaseList = cellstr(phaseList);
end
  

%% Determine type of pick file
if DFS(end) == '/'
  DFS = DFS(1:end-1);
end
if exist(DFS) == 7
  if exist([DFS '/tlControl.mat'])==2
    pickFormat = 3;
  else
    pickFormat = 1;
  end
elseif exist(DFS)==1
  pickFormat = 2;
elseif exist(DFS) || exist([DFS '.mat']);
  pickFormat = 2;
else
  status = 2;
  return
end

% Get picks
% tlPick ASCII files
if pickFormat == 1;
  if isempty(chanIter)
    channelSpecific = false;
  elseif isnan(chanIter)
    channelSpecific = false;
  elseif ~chanIter
    channelSpecific = false;
  elseif chanIter<0
    channelSpecific = true;
  elseif chanIter>0
    tMD.channel = zeros(size(tMD.channel))+chanIter;
    channelSpecific = true;
  end
%   if ~iscell(phaseList);
%     phaseList = string2cell(phaseList);
%   end
  if isempty(phaseList)
    status = 4;
    return
  else
    pick = [];
    for i = 1:length(phaseList)  
      [tlPick, status] = load_tlPick(DFS, tMD, phaseList{i}, channelSpecific, {'time'});
      if ~status
        pick = [pick tlPick.time];
      else
        pick = [pick nan(tMD.ntrace,1)];
      end
    end
  end
    
% tlArrival file 
elseif pickFormat == 2
  [pick, phase, status] = load_tlArrivalPickMatrix(DFS, tMD, phaseList);
  if status
    status = status + 1;
  end
  
% tlMisfit file
elseif pickFormat == 3
  if isempty(chanIter)
    chanIter = 1;
  elseif chanIter==0
    chanIter = 1;
  end
  [pick, phase, status] = load_tlOutputPickMatrix(DFS, tMD, chanIter, phaseList);
  if status
    status = status + 1
  end
end
    
  
  
