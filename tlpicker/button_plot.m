% Script to plot record section

menu_off
drawnow

if ~doneWhat.segy
  beep
  warning('Failure - No SEGY loaded to plot')

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
  
  % Verify that something has changed
  if upToDate.static && upToDate.mute && upToDate.sort && upToDate.inactive && upToDate.plot && doneWhat.plot && ~loadingView
    warning('The current plot is already up to date')

  else   
    if ~loadingView
      % Load static corrections
      if ~upToDate.static
        if isempty(deblank(currentMenu.static.DFS));
          traceMetaData.static = zeros(traceMetaData.ntrace,1);      
        else
          [static, status] = get_pickMatrix(currentMenu.static.DFS, ...
                                 currentMenu.static.phase, ...
                                 currentMenu.static.chanIter, traceMetaData);
          if ~status
            if size(static,2) ~= 1
              traceMetaData.static = zeros(traceMetaData.ntrace,1);      
              warning('Static options cannot define multiple picks per trace - No statics')
            else
              if any(isnan(static))
                static(isnan(static)) = 0;
                warning('some traces have no static correction - Set to 0 for these traces')
              end
              traceMetaData.static = static;
            end
          else
            if status==2
              traceMetaData.static = zeros(traceMetaData.ntrace,1);  
              warning(['Static DFS ' currentMenu.static.DFS ' is not a valid directory/file/structure - No statics'])
            elseif status==3
              traceMetaData.static = zeros(traceMetaData.ntrace,1);  
              warning(['Static DFS ' currentMenu.static.DFS ' contains multiple picks for a station/eventid/channel/phase - No statics'])
            elseif status==4
              traceMetaData.static = zeros(traceMetaData.ntrace,1);  
              warning('Static phase needed for a tlPick file - No statics')
            end   
          end
        end
        upToDate.static = true;
      end

      % Load mute times
      if ~upToDate.mute
        if isempty(deblank(currentMenu.mute.DFS))
          traceMetaData.mute = Inf(traceMetaData.ntrace,1);
        else
          [mute, status] = get_pickMatrix(currentMenu.mute.DFS, ...
                                 currentMenu.mute.phase, ...
                                 currentMenu.mute.chanIter, traceMetaData);      
          if ~status
            if size(mute,2) ~= 1
              traceMetaData.mute = Inf(traceMetaData.ntrace,1);
              warning('Mute options cannot define multiple picks per trace - No mute')
            else
              if any(isnan(mute))
                mute(isnan(mute)) = Inf;
                warning('some traces have no mute correction - Set to Inf for these traces')
              end
              traceMetaData.mute = mute;
            end
          else
            traceMetaData.mute = Inf(traceMetaData.ntrace,1); 
            if status==2
              warning(['Mute DFS ' currentMenu.mute.DFS ' is not a valid directory/file/structure - No mute'])
              traceMetaData.mute = Inf(traceMetaData.ntrace,1);
            elseif status==3
              warning(['Mute DFS ' currentMenu.mute.DFS ' contains multiple picks for a station/eventid/channel/phase - No mute'])
              traceMetaData.mute = Inf(traceMetaData.ntrace,1);
            elseif status==4
              warning('Mute phase needed for a tlPick file - No mute')
              traceMetaData.mute = Inf(traceMetaData.ntrace,1);
            end
          end
        end                    
        upToDate.mute = true;
      end

      % Load inactive picks
      if ~upToDate.inactive
        inactive.time = [];
        inactive.lineType = {};
        if isempty(deblank(currentMenu.inactive.DFS))

        else
          phaseList = string2cell(currentMenu.inactive.phase);
          for DFS = string2cell(currentMenu.inactive.DFS);
            [pick, status] = get_pickMatrix(DFS, phaseList, ...
                          currentMenu.inactive.chanIter, traceMetaData);
            if ~status
              inactive.time = [inactive.time pick];
            else
              if status==2
                warning(['Inactive DFS ' DFS ' is not a valid directory/file/structure'])
              elseif status==3
                warning(['Inactive DFS ' DFS ' contains multiple picks for a station/eventid/channel/phase - No mute'])
              elseif status==4
                warning('A phase is needed for inactive tlPick files')
              end
            end
          end
          lineType = string2cell(currentMenu.inactive.lineType);
          if isempty(lineType)
            lineType = {''};
          end
          for i = 1:size(inactive.time,2);
            inactive.lineType(i) = lineType(rem(i-1,length(lineType))+1);
          end
          upToDate.inactive = true;
        end
      end

      % Trace sorting
      if ~upToDate.sort
        sortIndex = get_sortIndex(traceMetaData,currentMenu.sort);
        traceMetaData.sortIndex = sortIndex;
        upToDate.sort = true;
      end
    end

    % Plotting here
    figure(hFig)
    clf
    [status,hRS,x0Trace,pickWidth,absScale] = plot_recordSection(traceData, traceMetaData, currentMenu);
    hold on
    lastPlot = currentMenu.plot;

    if ~status
      % Plot inactive picks
      if ~isempty(inactive.time);
        for i = 1:size(inactive.time,2);
          plot_pick(inactive.time(:,i), x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, inactive.lineType{i});
        end
      end
      doneWhat.plot = true;
      upToDate.plot = true;
      
      % Plot inactive picks of record section
      plot_geometryInactive
      
      % Plot active picks if they exist
      if doneWhat.loadPick
        figure(hFig)
%        [status, handlesActive.time] = plot_pick(active.time, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r',active.use);
        [status, handlesActive.time] = plot_pick(active.time, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r-2',active.use);
        [status, handlesActive.unc1] = plot_pick(active.time - active.unc, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r:',active.use);
        [status, handlesActive.unc2] = plot_pick(active.time + active.unc, x0Trace, traceMetaData.range, lastPlot.redVel, traceMetaData.static, pickWidth, 'r:',active.use);
        plot_geometryActive
      end
    elseif status == 1
      warning('Invalid X-axis Option - No plot ')
    elseif status == 2
      warning('No traces within X limits - No plot ')
    elseif status == 3
      warning('No data within time limits - No plot ')
    end
  end
end 

menu_on
drawnow
