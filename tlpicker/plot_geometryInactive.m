% Plots all shots with at least one inactive pick on experiment basemap

if length(unique(traceMetaData.channel))>1 || length(unique(traceMetaData.station))>1
  warning('Cannot plot events with inactive picks on experiment geometry plot when there is more than one station or channel')
else
  
  if exist('hGP3')==1
    figure(100)
    delete(hGP3)
    clear hGP3
  end
  
  if currentMenu.geometryPlot.inactive
  
    tMDtemp.ntrace = srEvent.nevt;
    tMDtemp.eventid = srEvent.id;
    tMDtemp.station = repmat(traceMetaData.station(1),tMDtemp.ntrace,1);
    tMDtemp.channel = repmat(traceMetaData.channel(1),tMDtemp.ntrace,1);

    iplot = false(tMDtemp.ntrace,1);
    phaseList = string2cell(currentMenu.inactive.phase);
    for DFS = string2cell(currentMenu.inactive.DFS);
       [temp, status] = get_pickMatrix(DFS, phaseList, ...
                          currentMenu.inactive.chanIter, tMDtemp);
       if ~status
         if size(temp,2)==1
           iplot = iplot | ~isnan(temp);
         else
           iplot = iplot | any(~isnan(temp'))';
         end
       end
    end

    if any(iplot)
      figure(100)
      hGP3 = plot(srEvent.x(iplot),srEvent.y(iplot),'ob','markersize',6);
      if exist('hGP1')==1
        legend('Event','Station','Loaded Event','Loaded Station','Inactive','location','bestoutside')
      else
        legend('Event','Station','Inactive','location','bestoutside')
      end
    end
  end
end
figure(hFig)
      