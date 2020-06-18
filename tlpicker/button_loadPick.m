% Script to load and plot active picks

% Turn off menu
menu_off
drawnow

% Need to have loaded SEGY
if ~doneWhat.segy
  beep
  warning('Failure - No SEGY loaded (need to load SEGY and plot a record section first')

% and to have plotted it
elseif ~doneWhat.plot
  beep
  warning('Failure - No record section plotted (need to plot record section first)')
  
else
  
  if ~upToDate.experiment
    warning('Experiment parameters have been modified since data loaded')
  end
  if ~upToDate.segy
    warning('SEGY parameters have been modified since data loaded')
  end
  if doneWhat.filter && ~upToDate.filter
    warning('Filter parameters have been modified since filter was applied')
  end
  if ~upToDate.static
    warning('Static parameters have been modified since data plotted')
  end
  if ~upToDate.mute
    warning('Mute parameters have been modified since data plotted')
  end
  if ~upToDate.sort
    warning('Sorting parameters have been modified since data plotted')
  end
  if ~upToDate.plot
    warning('Plotting parameters have been modified since data plotted')
  end
  if ~upToDate.inactive
    warning('Inactive pick parameters have been modified since data plotted')
  end
  if doneWhat.madePick
    % Prompt to save picks if there are some
    warning('Loading picks will discard existing picks')
    proceed = input('Do you want to save existing picks first (y/N) : ','s');
    if strcmpi(proceed,'y')
      button_savePick
    end
    clear active
    figure(hFig)
    if exist('handlesActive')
      delete(handlesActive.time(isgraphics(handlesActive.time')));
      delete(handlesActive.unc1(isgraphics(handlesActive.unc1')));
      delete(handlesActive.unc2(isgraphics(handlesActive.unc2')));
      clear handlesActive
    end
  end
  
  [active, status] = load_tlPick(currentMenu.active.directory, ...
                       traceMetaData, currentMenu.active.phase , ...
                       currentMenu.active.channelSpecific, {'all'});
  
  if status==1
    warning('Failure - tlPick directory does not exist')
    doneWhat.loadPick = false;
  
  else
    if ~any(active.picked)
      warning('Warning: tlPick directory contains no picks for current traces and picks')
      disp('         Picking will proceed from scratch')
    end
    
    % Plot the picks and uncertainties (currently status is always false)
    figure(hFig)
%   [status, handlesActive.time] = plot_pick(active.time, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r',active.use);
    [status, handlesActive.time] = plot_pick(active.time, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r-2',active.use);
    [status, handlesActive.unc1] = plot_pick(active.time - active.unc, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r:',active.use);
    [status, handlesActive.unc2] = plot_pick(active.time + active.unc, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r:',active.use);
    plot_geometryActive
    doneWhat.loadPick = true;
    doneWhat.madePick = false;
    
  end
  
end

menu_on
drawnow    
  
