% Script to load SEGY data (uses catalogs of trace offsets)

% Turn off the menu
menu_off 
drawnow

% Update experiment
if ~upToDate.experiment
  srGeometry = load_srGeometry(currentMenu.experiment.srGeometry);
  srStation = load_srStation(currentMenu.experiment.srStation, srGeometry);
  srEvent = load_srEvent(currentMenu.experiment.srEvent, srGeometry);
  eval(['load ' currentMenu.experiment.segyCatalog]);    
  srStationRotation = load_srStationRotation(currentMenu.orientHoriz.srStationRotation,srStation);
  upToDate.experiment = true;
  mustLoad = true;
else
  mustLoad = false;
end

% Load data if necessary
if mustLoad || ~doneWhat.segy || ~upToDate.segy
  
  % Prompt to save picks if there are some
  if doneWhat.madePick
    warning('Loading SEGY will discard existing picks')
    proceed = input('Do you want to save existing picks first (Y/n) : ','s');
    if ~strcmpi(proceed,'n')
      button_savePick
    end
    clear active
  end
  
  
  [traceMetaData,traceData] = ...
         get_segy(currentMenu.segy, segyCatalog, ...
                  srEvent, srStation, srGeometry, srStationRotation);


  if isempty(traceData)
    warning('Failure - Load SEGY returned no data')
    doneWhat.segy = false;
    upToDate.segy = false;

  else
    fprintf('Loaded Data for Events\n')
    fprintf('%i %i %i %i %i %i %i %i %i %i %i %i\n'  , sort(unique(traceMetaData.eventid)))
    fprintf('\n')

    traceDataUnfilt = traceData; 
    constantSamprate = ~all(diff(traceMetaData.samprate)) | traceMetaData.ntrace==1;
    doneWhat.segy = true;
    upToDate.segy = true;
    upToDate.mute = false;
    upToDate.static = false;
    upToDate.sort = false;
    upToDate.inactive = false;
    upToDate.plot = false;
    
    % Plot Experiment Geometry
    plot_geometryLoaded

  end

  % Status of downstream operations
  doneWhat.filter = false;
  doneWhat.plot = false;
  doneWhat.loadPick = false;
  clear handlesActive
  doneWhat.madePick = false;
  figure(hFig)
  clf
  
  % Filter if done on load
  if currentMenu.filter.onLoad
    button_filterSEGY
  end
  
else
  warning('SEGY has already been loaded')
  
end

menu_on
drawnow
