% Script to verify changes in (active) picking menu for fields that are normally fixed
% Executed by callbacks in tlPicker_menu

menu_off
drawnow

if ~strcmp(fixedActive.directory,currentMenu.active.directory) && doneWhat.loadPick
  warning('The picking directory should not normally be changed during a picking session.')
  disp('Proceeding will lead to the status of all active picks being changed')
  disp('to updated but picks will not be automatically reloaded from the new directory')
  proceed = input('Do you wish to proceed (y/N) : ','s');
  if strcmpi(proceed,'y')
    active.updated(:) = true;
    fixedActive.directory = currentMenu.active.directory;
  else
    disp('Picking directory not changed')
    currentMenu.active.directory = fixedActive.directory;
    set(handlesMenu.active_directory,'string',currentMenu.active.directory)
  end
end

if ~strcmp(fixedActive.phase,currentMenu.active.phase) && doneWhat.loadPick
  if doneWhat.madePick
    warning('The picking phase cannot be changed when picks are loaded')
    proceed = input('Do you wish to proceed and discard current active picks without saving (y/N) : ','s');
  else
    warning('The picking phase cannot be changed when picks are loaded')
    proceed = input('Do you wish to proceed and discard current active picks (y/N) : ','s');
  end
  if strcmpi(proceed,'y')
    doneWhat.loadPick = false;
    doneWhat.madePick = false;
    active = [];
    delete(handlesActive.time(isgraphics(handlesActive.time)))
    delete(handlesActive.unc1(isgraphics(handlesActive.unc1)))
    delete(handlesActive.unc2(isgraphics(handlesActive.unc2))) 
    clear handlesActive
    fixedActive.phase = currentMenu.active.phase;
    disp('You can now load the new picks');
  else
    disp('Picking phase not changed')
    currentMenu.active.phase = fixedActive.phase;
    set(handlesMenu.active_phase,'string',currentMenu.active.phase)
  end
end

if fixedActive.channelSpecific~=currentMenu.active.channelSpecific && doneWhat.loadPick
  warning('Changing the Channel Specific designation during pickingis unsupported and may have bad consequences')
  proceed = input('Do you wish to proceed and change the "Channel Specific" (y/N) : ','s');
  if strcmpi(proceed,'y')
    fixedActive.channelSpecific = currentMenu.active.channelSpecific;
  else
    disp('Picking "Channel Specific" option has not changed')
    currentMenu.active.channelSpecific = fixedActive.channelSpecific;
    set(handlesMenu.active_channelSpecific,'value',currentMenu.active.channelSpecific)
  end
end

menu_on
drawnow

