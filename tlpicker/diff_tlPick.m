function tdiff = diff_tlPick(srEvent, station, dirName1, onephase1, dirName2, onephase2)
% Returns difference in picks for tlPick files for a single/multiple stations
%
% Usage
%   tdiff = diff_tlPick(srEvent, station, dirName1, phase1, dirName2, phase2)
%
% Inputs
%   srEvent     - Stingray Event structure (see Stingray manual)
%   station     - Station name(s) in cell structure
%   dirName1    - Directory of 1st tlpick file
%   onephase1   - Phase name of 1st tlpick file
%   dirName2    - Directory of 2nd tlpick file
%   onephase2   - Phase name of 2nd tlpick file
%
% Outputs
%   tdiff       - Matrix of pick times (pick1 - pick2) with NaN for no data
%                 and one column for each station
%
% Presently this function will not deal with channel specific picks

% Process inputs
if ischar(station)
  station = {station};
end

tdiff = NaN(srEvent.nevt,length(station));

for ista = 1:length(station)
  fprintf(['diff_tlPick processing station ' station{ista} '\n']);

  [station1, ~, eventid1, phase1, time1] = read_tlPick(dirName1, station{ista}, onephase1);
  [station2, ~, eventid2, phase2, time2] = read_tlPick(dirName2, station{ista}, onephase2);
  
  if ~isempty(station1) && ~isempty(station2)

      if any(~strcmp(station1, station{ista}))
        error('First tlPick file has unexpected station')
    %     tdiff = [];
    %     return
      end
      if any(~strcmp(station2, station{ista}))
        error('Second tlPick file has unexpected station')
    %     tdiff = [];
    %     return
      end
      if any(~strcmpi(phase1, onephase1))
        error('First tlPick file has unexpected phase')
    %     tdiff = [];
    %     return
      end
      if any(~strcmpi(phase2, onephase2))
        error('Second tlPick file has unexpected phase')
    %     tdiff = [];
    %     return
      end

      for i = 1:length(eventid1);
        j = find(eventid2 == eventid1(i));
        if length(j)>1
          error('Repeated picks in second tlPick file')
        elseif ~isempty(j)
          k = find(srEvent.id == eventid1(i));
          if isempty(k)
            error('Event ID not in srEvent')
          else
            tdiff(k,ista) = time1(i) - time2(j);
          end
        end
      end
  end
end



