
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
    end