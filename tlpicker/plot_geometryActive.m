% Plots all shots with an active pick on experiment basemap

if length(unique(traceMetaData.channel))>1 || length(unique(traceMetaData.station))>1
  warning('Cannot plot events with active picks on experiment geometry plot when there is more than one station or channel')
else
 
  if exist('hGP4')==1
    figure(100)
    delete(hGP4)
    clear hGP4
  end
  
  if currentMenu.geometryPlot.active

    tMDtemp.ntrace = srEvent.nevt;
    tMDtemp.eventid = srEvent.id;
    tMDtemp.station = repmat(traceMetaData.station(1),tMDtemp.ntrace,1);
    tMDtemp.channel = repmat(traceMetaData.channel(1),tMDtemp.ntrace,1);


    [temp, status] = get_pickMatrix(currentMenu.active.directory, currentMenu.active.phase, ...
                      currentMenu.active.channelSpecific, tMDtemp);
    if ~status
      if size(temp,2)==1
        iplot = ~isnan(temp);
      else
        iplot = any(~isnan(temp'))';
      end
    end

    if any(iplot)
      figure(100)
      hGP4 = plot(srEvent.x(iplot),srEvent.y(iplot),'og','markersize',8);
      if exist('hGP1')==1 && exist('hGP3')==1
        legend('Event','Station','Loaded Event','Loaded Station','Inactive','Active','location','bestoutside')
      elseif exist('hGP1')==1
        legend('Event','Station','Loaded Event','Loaded Station','Active','location','bestoutside')
      elseif exist('hGP3')==1
        legend('Event','Station','Inactive','Active','location','bestoutside')
      else
        legend('Event','Station','Active','location','bestoutside')      
      end
     
    end
  end
end  
figure(hFig)