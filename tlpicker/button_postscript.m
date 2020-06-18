% Button to save the record section to a Postscript file

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
  warning('Failure - No record section plotted (need plot record section first')
  
else
  if isempty(currentMenu.ps.fileName)
    if currentMenu.ps.append
      fileName = ['tlPicker_plot_' date '.eps'];
    else
      fileName = ' ';
      i = 0;
      while exist(fileName)==2 || strcmp(fileName,' ');
        i = i + 1;
        fileName = ['tlPicker_plot_' date '_' int2str(i) '.eps'];
      end
    end
  else
    fileName = currentMenu.ps.fileName;
  end
  if length(fileName)<5
    fileName = [fileName '.eps'];
  end
  if ~strcmp(fileName(end-3:end),'.eps')
    fileName = [fileName '.eps'];
  end
  if currentMenu.ps.append
    drawnow
    eval(['print -f' int2str(hFig) ' -dpsc2 -append ' fileName]);
  else
    if exist(fileName)==2
      warning(['Warning: Postscript file ' fileName ' will be overwritten']);
    end
    drawnow
    eval(['print -f' int2str(hFig) ' -dpsc2 ' fileName]);
  end
  disp(['Current record section plot saved to file ' fileName]);
  
end

menu_on
drawnow