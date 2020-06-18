% Script to load the current view

menu_off
drawnow

if isempty(currentMenu.view.fileName)
  beep
  warning('Failure: no view will be loaded because the file name is not specified')

else
  if doneWhat.madePick
    beep
    warning(' There are currently unsaved picks')
    disp('Unsaved picks will be lost by loading a view')
    proceed = input('Do you want to proceed (y/N) : ','s');
    if strcmpi(proceed,'y')
      proceed = true;
    else
      proceed = false;
    end
  else
    proceed = true;
  end
  
  if proceed   
    if ~exist(currentMenu.view.fileName)==2 || ...
       ~exist([currentMenu.view.fileName '.mat'])==2
      beep
      warning('Failure: no view will be loaded because the file does not exist')
    
    else
      
      eval(['load ' currentMenu.view.fileName])
      currentMenu.experiment = defaultMenu.experiment;
      tempDoneWhat = doneWhat;
      clear doneWhat
      tlPicker
      if tempDoneWhat.segy
        doneWhat.segy = false;
        button_loadSEGY
        traceMetaData.static = static;
        traceMetaData.mute = mute;
        traceMetaData.sortIndex = sortIndex;
        doneWhat = tempDoneWhat;
        if doneWhat.plot
          loadingView = true;
          button_plot
          loadingView = false;
        end
      end
      
    end
    
  end
  
end
    
menu_on
drawnow    