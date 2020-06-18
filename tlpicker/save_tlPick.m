function status = save_tlPick(dirName, tlPick, channelSpecific, orderEventid)
% Save a tlPick structure to ASCII tlPick_STATION_PHASE.dat files
%
% Usage
%   status = save_tlPick(dirName, tlPick, channelSpecific, orderEventid)
%
% Inputs
%   dirName         - Directory for tlPick ASCII
%   tlPick          - tlPick structure
%                     Required fields that are saved are
%                       station
%                       channel
%                       eventid
%                       phase
%                       time
%                       unc
%                     Optional fields that are also saved are
%                       filtLim0
%                       filtLim1
%                       filtOrder
%                       filtZeroPhase
%                       scale
%                       user
%                       lddate
%                       use
%                       comment
%                     Optional fields that are not saved are
%                       stationList
%                       phaseList - Always only 1 phase for tlPicker
%                       picked
%                       updated          
%   channelSpecific - Logical to control how existing picks are replaced
%                       FALSE - Channel is not considered
%                               A pick is replaced if it has a 
%                               station/eventid/phase match
%                               An error will occur if an ASCII file has 
%                               more than one pick for a station/eventid/phase 
%                               match (even with different channels)
%BUG - At present no checkin of tlPICK to ensure that it has no duplicates
%                        TRUE - Channel is considered
%                               A pick is replaced if the 
%                               station/event/phase/channel matches
%   orderEventid    - Logical to control how picks are replaced or deleted
%                       FALSE - When picks are deleted the pick phase is
%                               marked as deleted and any new pick appended
%                               to the end of the ASCII file
%                       TRUE  - The whole file is rewritten if necessary to
%                               eliminate deleted picks and order picks by
%                               eventid and then channel.
%
% Outputs       
%   status - Status of execution
%            0 - Okay
%            1 - Directory does not exist
%            2 - One of the required tlPick fields does not exist

global debugPicking

status = 0;
lengthRecord = 213; % Length of record in tlPick file

%% Process inputs
if nargin<3
  channelSpecific = 0;
end

if nargin<4
  orderEventid = 0;
end

if ~isempty(dirName)
  if dirName(end) == '/'
    dirName = dirName(1:end-1);
  end
else
  dirName = '.';
end
if exist(dirName)~=7
  warning(['save_tlPick: Directory ' dirName ' does not exist']);
  status = 1;
  return
end

if ~all(isfield(tlPick,{'station','channel','eventid','phase','time','unc'}))
  warning(['save_tlPick: At least one required tlPick fields does not exist']);
  status = 2;
  return
end

%% Create pick variables from pick structure (speeds things up)
npick = length(tlPick.time);
station = tlPick.station;
channel = tlPick.channel;
eventid = tlPick.eventid;
phase = tlPick.phase;
time = tlPick.time;
unc = tlPick.unc;
if ~isfield(tlPick,'filtLim0')
  filtLim0 = zeros(npick,1);
else
  filtLim0 = tlPick.filtLim0;
end
if ~isfield(tlPick,'filtLim1')
  filtLim1 = zeros(npick,1);
else
  filtLim1 = tlPick.filtLim1;
end
if ~isfield(tlPick,'filtOrder')
  filtOrder = zeros(npick,1);
else
  filtOrder = tlPick.filtOrder;
end
if ~isfield(tlPick,'filtZeroPhase')
  filtZeroPhase = zeros(npick,1);
else
  filtZeroPhase = tlPick.filtZeroPhase;
end
if ~isfield(tlPick,'scale')
  scale = zeros(npick,1);
else
  scale = tlPick.scale;
end
if ~isfield(tlPick,'user')
  user = repmat({'unknown'},npick,1);
else
  user = tlPick.user;
end
if ~isfield(tlPick,'lddate')
  c = clock;
  lddate = repmat(date2secnds(c(1),julday(c(2),c(3),c(1)),c(4),c(5),c(6)),npick,1);
else
  lddate = tlPick.lddate;
end
if ~isfield(tlPick,'use')
  use = true(npick,1);
else
  use = tlPick.use;
end
if ~isfield(tlPick,'comment')
  comment = repmat({'None'},npick,1);
else
  comment = tlPick.comment;
end

% Lists of stations and picks in file
if ~isfield(tlPick,'stationList')
  stationList = unique(tlPick.station);
else
  stationList = tlPick.stationList;
end
if ~isfield(tlPick,'phaseList')
  phaseList = unique(tlPick.phase);
else
  phaseList = tlPick.phaseList;  
end

% Default status of picks
if ~isfield(tlPick,'picked')
  picked = ~isnan(time);
else
  picked = tlPick.picked;
end
if ~isfield(tlPick,'updated')
  updated = true(npick,1);
else
  updated = tlPick.updated;
end

% Loop over all stations and phases 
for currentStation = stationList

  for currentPhase = phaseList
  
    % Test for picks for station / phase combination
    iout = find(strcmp(tlPick.station,currentStation) & strcmp(tlPick.phase,currentPhase));
    if ~isempty(iout)
      
      % tlPick file name
      fileName = [dirName '/tlPick_' currentStation{1} '_' currentPhase{1} '.dat'];

      %% Process content of existing file
      if exist(fileName)
        
       % Get their details
        [savedStation,savedChannel,savedEventid,savedPhase,savedTime] = ...
          read_tlPick(dirName, currentStation{1}, currentPhase, ...
                      {'station','channel','eventid','phase','time'});

 
        if ~isempty(savedStation)
          
          % Verify station and phases are correct
          if ~all(strcmp(savedStation,currentStation)) && ...
             ~all(strcmp(savedPhase,currentPhase))      
            warning(['save_tlPick: file ' fileName ' has mismatched station or phase'])
            error('Aborting')
          end
        
          % Logicals to indicate whether picks are to be kept or not
          keep = true(length(savedEventid),1);
          if orderEventid
            keep(isnan(savedTime)) = false;
          end

          % Loop through new picks
          for i=1:length(iout);

            if updated(iout(i))

              % Find index of old undeleted pick corresponding to new picks or deletes
              if channelSpecific
                j = find(channel(iout(i))==savedChannel & eventid(iout(i))==savedEventid);
              else
                j = find(eventid(iout(i))==savedEventid);
              end
              j = j(~isnan(savedTime(j)));

              % Tag old pick to be deleted
              if length(j)==1
                keep(j) = false;
              elseif length(j)>1
                warning(['save_tlPick: file ' fileName ' has duplicate picks'])
                error('Aborting')
              end

             end

          end
          
        end
      
      end
      
      %% Deal with existing pick file and then open it ready to write
      if ~exist(fileName) || isempty(savedStation)       
      % No existing file or it is empty
        c1 = [];
        eventidSort = [];
        channelSort = [];
        fid = fopen(fileName,'w');
        
      elseif orderEventid && any(keep)        
        % Existing data be resorted  
        fid = fopen(fileName,'r');
        [c1,n] = fread(fid,inf,'uint8');
        fclose(fid);
        c1 = reshape(c1,lengthRecord,n/lengthRecord);
        c1 = c1(:,keep);
        eventidSort = savedEventid(keep);
        channelSort = savedChannel(keep);
        fid = fopen(fileName,'w');

      elseif orderEventid
      % Resorting but no existing data to keep 
        c1 = [];
        eventidSort = [];
        channelSort = [];      
        fid = fopen(fileName,'w');
        
      elseif ~orderEventid & all(keep)
      % No resorting and no deleting of existing data
        c1 = [];
        fid = fopen(fileName,'a');
      
      % No resorting and deleting of existing data
      elseif ~orderEventid
        c1 = [];
        fid = fopen(fileName,'r+');
        for i = find(~keep(:)');
          fseek(fid,(i-1)*lengthRecord+16,'bof')
          fwrite(fid,abs(' Deleted'),'uint8')
        end
        fseek(fid,0,'eof')
          
      end 
      
      
      %% Create new pick output     
      c2 = zeros(1,lengthRecord*sum(picked(iout) & updated(iout)),'uint8');
      j = 0;
      for i = iout(:)'
        if debugPicking  && ~channelSpecific
          if any(eventid(i)==eventid([1:i-1 i+1:end]))
            pause(0.5); beep; pause(0.5); beep; pause(0.5); beep; pause(0.5); beep; pause(0.5); beep; pause(0.5); beep 
            disp('Potential bug with picking')
            disp('Workspaces have been saved in matlab_main4William.mat and matlab_save_tlPick4William')
            disp('Send these files and the pick file to William')
            disp('Type return to continue')
            evalin('base', 'menu_off')
            evalin('base', 'save matlab_main4William')
            save matlab_save_tlPick4William
            keyboard
            evalin('base', 'menu_on')
          end
        end
        if picked(i) && updated(i);
          c2(j+1:j+lengthRecord) = ...
            sprintf('%6s%2i%8i%8s%10.4f%7.4f%5.1f%5.1f%2i%2i%10.2e%10s%17.5f%2i%38s%-80s\n', ...
                  station{i}, channel(i), eventid(i), phase{i}, ...
                  time(i), unc(i), filtLim0(i), filtLim1(i), ...
                  filtOrder(i), filtZeroPhase(i), scale(i), user{i}, ...
                  lddate(i),use(i), ' ', comment{i});
          j = j + lengthRecord;
        end
      end
      
      %% Order output
      if orderEventid
        eventidSort = [eventidSort; eventid(iout(picked(iout) & updated(iout)))];
        channelSort = [channelSort; channel(iout(picked(iout) & updated(iout)))];
        c2 = [c1 reshape(c2,lengthRecord,sum(picked(iout) & updated(iout)))];
        [channelSort,i] = sort(channelSort);
        eventidSort = eventidSort(i);
        c2 = c2(:,i);
        [eventidSort,i] = sort(eventidSort);
        c2 = c2(:,i);
      end

      
      %% Write Data 
      fwrite(fid, c2, 'uint8');
      fclose(fid);
      
    end % if any(iout)
    
  end % currentPhase
  
end  % currentStation
 

  

  