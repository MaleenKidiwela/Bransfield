% Button to print the record section

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
  figure(hFig)
  orient landscape
  print
end

menu_on
drawnow