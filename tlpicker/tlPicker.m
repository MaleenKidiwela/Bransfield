% Script to run TomoLab Picker
% 2/9/2020 - Need to turn off toolbar in figures to rubberbanding works

% Added to try and understand an occasional error with duplicate picks by
% saving graphical inputs and also checking pick structure
global debugPicking
debugPicking = false;
if debugPicking
  ginputRecord = [];
end
  
% Do not allow tlPicker to be started from within another function (except
% for button_loadView
temp = dbstack;
if ~(length(temp)==1 || (length(temp)==2 && strcmp(temp(2).name,'button_loadView')) ...
    || (length(temp)==2 && strcmp(temp(2).name,'fix_sections')))
  error('tlPicker has to be started from the main workspace (i.e., no ''K>>'' prompt)')
end

% Create default menu values unless restarting from a crash
if exist('currentMenu')~=1
  if exist('default_menu.m')==2
    default_menu
  else
    default_menu_example
  end
  currentMenu = defaultMenu; 
end

% Save some fields that should not normally be changed
fixedActive.directory = currentMenu.active.directory;
fixedActive.phase = currentMenu.active.phase;
fixedActive.channelSpecific = currentMenu.active.channelSpecific;

% Clear figures (100 - Experiment geometry, Starting 101 Figure Window)
if ~exist('hFig','var'); hFig=101; end 
figure(100); clf;
set(100,'toolbar','none')
figure(hFig); clf
set(hFig,'toolbar','none')
clear hGP1 hGP2 hGP3 hGP4

% Turn of view loading flag for button_plot
loadingView = false;

% Set status structures
if ~exist('doneWhat')
  doneWhat.segy = false;
  doneWhat.filter = false;
  doneWhat.plot = false;
  doneWhat.loadPick = false;
  doneWhat.madePick = false;
  upToDate.experiment = false;
  upToDate.segy = false;
  upToDate.filter = false;
  upToDate.mute = false;
  upToDate.static = false;
  upToDate.sort = false;
  upToDate.inactive = false;
  upToDate.plot = false;  
  upToDate.picking = false;
else
  % Restarting obliterates record section
  doneWhat.plot = false;
  if doneWhat.segy
    plot_geometryLoaded
  end
end

% Create the menu with gui handles in the main workspace
handlesMenu = tlPicker_menu_export(currentMenu);
% handlesMenu = tlPicker_menu(currentMenu);
menu_on

% Disable the keyboard
% uiwait(handles.figure1)
