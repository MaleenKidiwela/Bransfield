% Script to save picks

if ~doneWhat.madePick
  if doneWhat.loadPick
    beep
    warning('Failure - There are no new picks to save')
  else
    beep
    warning('Failure - There are no picks to save')
  end
else
  if exist(currentMenu.active.directory)~=7
    warning(['Directory ' currentMenu.active.directory ' does not exist'])
    proceed = input('Do you want to create it (y/N) :','s');
    if strcmpi(proceed,'y') 
      proceed = true;
      unix(['mkdir ' currentMenu.active.directory])
    else
      proceed = false;
    end
  else 
    proceed = true;
  end
  if proceed
    status = save_tlPick(currentMenu.active.directory, active, ...
                         currentMenu.active.channelSpecific, ...
                         currentMenu.active.orderEventID);
    if status
      warning('Something went wrong with saving tlPick')
      keyboard
    else
      doneWhat.madePick = false;
      active.updated(:) = false;
      % Check to make sure that saved picks are not used as static
      % corrections, mute times or inactive picks (if so they will need to
      % be reloaded when data is plotted).  There is no such check for
      % currentMenu.segy.pickSelect
      if any(strcmp(string2cell(currentMenu.inactive.DFS),currentMenu.active.directory)) && ...
         any(strcmp(string2cell(currentMenu.inactive.phase),currentMenu.active.phase))
        upToDate.inactive = false;
      end
      if strcmp(currentMenu.static.DFS,currentMenu.active.directory) && ...
         strcmp(currentMenu.static.phase,currentMenu.active.phase)
        upToDate.static = false;
      end
      if strcmp(currentMenu.mute.DFS,currentMenu.active.directory) && ...
         strcmp(currentMenu.mute.phase,currentMenu.active.phase)
        upToDate.mute = false;
      end
      plot_geometryActive
      disp('Picks have been saved')          
      
    end
  end 
end
      
                       
