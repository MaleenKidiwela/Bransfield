% Script to filter (or refilter) SEGY data that has already been loaded

menu_off
drawnow

if ~doneWhat.segy
  beep
  warning('Failure - No SEGY loaded to filter')

else
  
  % Only warnings if upstream parameters have been changed
  if ~upToDate.experiment
    warning('Experiment parameters have been changed')
  end
  if ~upToDate.segy
    warning('SEGY parameters have been changed')
  end
  
  if ~doneWhat.filter || ~upToDate.filter

    success = 1;
    flim = zeros(1,2);
    if isempty(currentMenu.filter.lim0)
      flim(1) = 0;
    elseif isnan(currentMenu.filter.lim0)
      flim(1) = 0;
    else      
      flim(1) = currentMenu.filter.lim0;
    end
    if isempty(currentMenu.filter.lim1)
      flim(2) = 0;
    elseif isnan(currentMenu.filter.lim1)
      flim(2) = 0;
    else
      flim(2) = currentMenu.filter.lim1;
    end
    if isempty(currentMenu.filter.order)
      warning('button_filterSEGY.m - no filter applied because order is undefined')
      success = 0;
    elseif isnan(currentMenu.filter.order)
      warning('button_filterSEGY.m - no filter applied because order is undefined')
      success = 0;
    elseif currentMenu.filter.order == 0
      warning('button_filterSEGY.m - no filter applied because order is zero')
      success = 0;
    else
      if constantSamprate
        wlim = flim * 2 * (1/traceMetaData.samprate(1));
      else
        wlim = (2 ./traceMetaData.samprate) * flim;
      end

      if any(any(wlim < 0))
        warning('button_filterSEGY.m - no filter applied because Fmin too low')
        success = 0;

      elseif any(any(wlim >= 1))
        warning('button_filterSEGY.m - no filter applied because Fmax too high')
        success = 0;

      elseif constantSamprate

        if wlim(1) && ~wlim(2)
          if currentMenu.filter.zeroPhase
            [b,a] = butter(ceil(currentMenu.filter.order/2),wlim(1),'high');
            traceData = filtfilt(b,a,traceDataUnfilt);
          else
            [b,a] = butter(ceil(currentMenu.filter.order),wlim(1),'high');
            traceData = filter(b,a,traceDataUnfilt);
          end 
          disp(['Filtered with a high-pass filter at ' num2str(flim(1)) ' Hz ']);

        elseif ~wlim(1) && wlim(2);
          if currentMenu.filter.zeroPhase
            [b,a] = butter(ceil(currentMenu.filter.order/2),wlim(2),'low');
            traceData = filtfilt(b,a,traceDataUnfilt);
          else
            [b,a] = butter(ceil(currentMenu.filter.order),wlim(2),'low');
            traceData = filter(b,a,traceDataUnfilt);
          end 
          disp(['Filtered with a low-pass filter at ' num2str(flim(2)) ' Hz '])

        elseif wlim(1) && wlim(2);
          if currentMenu.filter.zeroPhase
            [b,a] = butter(ceil(currentMenu.filter.order/2),wlim);
            traceData = filtfilt(b,a,traceDataUnfilt);
          else
            [b,a] = butter(ceil(currentMenu.filter.order),wlim);
            traceData = filter(b,a,traceDataUnfilt);
          end 
          disp(['Filtered with a band-pass filter at [' num2str(flim) '] Hz '])

        else
          warning('No filter applied because neither limit is set')

        end

      else
        warning('button_filterSEGY.m - variable sampling rate filtering not implemented')
        success = 0;
      end

      % Mute non zero phase
      if success && ~currentMenu.filter.zeroPhase
        if ~isempty(currentMenu.filter.muteTime)
          if currentMenu.filter.muteTime>0
            nmute = round(traceMetaData.samprate * currentMenu.filter.muteTime);
            for i = 1:traceMetaData.ntrace
              traceData(1:min(traceMetaData.nsamp(i),nmute(i)),i) = NaN;
            end
          end
        end
      end
    end

    if ~success
      doneWhat.filter = false;
      upToDate.filter = false;

    else
      doneWhat.filter = true;
      upToDate.filter = true;
      lastFilter = currentMenu.filter;
      if isempty(lastFilter.lim0)
        lastFilter.lim0 = 0;
      elseif isnan(lastFilter.lim0)
        lastFilter.lim0 = 0;
      end
      if isempty(lastFilter.lim1)
        lastFilter.lim1 = 0;
      elseif isnan(lastFilter.lim1)
        lastFilter.lim1 = 0;
      end
      if isempty(lastFilter.order)
        lastFilter.order = 0;
      elseif isnan(lastFilter.order)
        lastFilter.order = 0;
      end

      if doneWhat.plot
        % If a plot exists we need to recreate it 
        % (this recursively plots picks)
        doneWhat.plot = false;
        button_plot
      end
        
      
    end

  else
    warning('Filtering has already been applied')

  end
end

menu_on
drawnow
    
  
