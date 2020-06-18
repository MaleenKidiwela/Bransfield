% Script to save the current view

menu_off
drawnow

if isempty(currentMenu.view.fileName)
  fileName = ' ';
  i = 0;
  while exist(fileName)==2 || strcmp(fileName,' ');
    i = i + 1;
    fileName = ['tlPicker_view_' date '_' int2str(i) '.mat'];
  end
else
  fileName = currentMenu.view.fileName;
end

if doneWhat.madePick
  warning('There are currently unsaved picks')
  disp('The view will only restore these if they are subsequently saved')
end

static = traceMetaData.static;
mute = traceMetaData.mute;
sortIndex = traceMetaData.sortIndex;
try
  active = active;
catch 
  active = [];
end
try
  inactive = inactive;
catch 
  inactive = [];
end
    
eval(['save ' fileName ' currentMenu doneWhat static mute sortIndex active inactive'])
disp(['Current view saved to file ' fileName]);
    
menu_on
drawnow    