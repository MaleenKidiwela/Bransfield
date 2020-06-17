function plot_obsloc_event_residual(s,ista,doMap,doNorm,lin_num,line_number)
% Preliminary function to create event residual plot(s) for obsloc inversions
% 
% Usage
%   plot_obsloc_event_residual(s,ista,doMap,doNormRms)
% Inputs
%   s         - obsloc output structure
%   ista      - Indicies of stations in s structure to plot 
%                (not station names or indicies in srStation input to obsloc)
%   doMap     - Logical 
%               0 - Make one plot of residuals versus event index
%                  color coded by station
%               1 - Make one map per station of color coded residuals
%               2 - Make one map for all stations of color coded residuals
%   doNorm    - Logical 0/1 to plot residual or residual normalized to
%               picking error
%   lin_num   - Logical 0/1 to plot for all lines or for one specific line
%   
%
% Could do enhancements such as offsetting of overlapping events and 
% labeling of events 


misfit = s.time - s.timePred;
%misfit = misfit(find((-0.04<misfit& misfit<0.04)));
if doNorm
  misfit = misfit ./ s.timeError;
end

if line_number
    line_number = line_number *1000;
end

if lin_num == 1
    for k = ista(1):ista(length(ista))
    find(~isnan(misfit(:,ista)));
    %staind = find(char(s.tlArrival.station) == char(int2str(ista)));
    eventind = find( s.srEvent.id>=line_number &  s.srEvent.id<(line_number + 1000));
    end
end
    
% Create a plot of station residuals versus station id
if ~doMap 
  %figure; clf
  style1 = ['bgrcmk']'; style2=['ox+*sdv^<>p']';
  style = [repmat(style1,11,1) repmat(style2,6,1)];
  nstyle = size(style,1);
  for j=1:length(ista)
    plot(s.ievt,misfit(:,ista(j)),style(rem(j-1,nstyle)+1,:)); hold on
  end
  if doNorm
    title('Normalized travel time misfits (Observed minus Predicted)','FontSize',11);
  else
    title('Travel time misfits (Observed minus Predicted)','FontSize',11);
  end
  xlabel('Event Index (obsloc input)'); 
  ylabel('Misfit, s');
  xlim([0 max(s.ievt)+1]);
  yplot = median(abs(misfit(~isnan(misfit))))*6;
  for j = 1:length(s.ievt)
    text(s.ievt(j),yplot,int2str(s.srEvent.id(j)), ...
             'verticalalignment','middle','rotation',90, ...
             'horizontalalignment','left','fontsize',8);
  end
  legend(s.srStation.name(ista),'location','best')
  
% Creates maps of residual
elseif doMap == 1
  for j=1:length(ista)
    figure
    if lin_num ==1
        data = misfit(eventind,ista(j))
    else
    data = misfit(:,ista(j));
    end
    lim = [min(data) max(data)];
    col = LinColor(data);
    plot(s.xEvent,s.yEvent,'.k','Markersize',5)
    hold on
    for i=1:length(data)      
      if ~isnan(data(i))
        plot(s.xEvent(i),s.yEvent(i),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8)
      end
    end
    plot(s.xStation(ista(j)),s.yStation(ista(j)),'^k','MarkerFaceColor','k')
    axis('equal');
    ylabel('Y, km');
    xlabel('X, km');
    title(['Station ' cell2mat(s.srStation.name(ista(j)))],'FontSize',11);
    caxis(lim)
    hcb = colorbar;
    h = gca;
    set(gcf,'currentaxes',hcb);
    if doNorm
      title('Norm Resid');
    else
      title('Residual, s');
    end
    set(gcf,'currentaxes',h);
  end

else
    % Plot all residuals on one map
    figure
    if lin_num ==1
        data = misfit(eventind,ista)
    else
        data = misfit(:,ista);
    end
    [m,n]=size(data); 
    data = reshape(data,m*n,1);
    lim = [min(data) max(data)];
    col = LinColor(data);
    plot(s.xEvent,s.yEvent,'.k','Markersize',5)
    hold on
    for i=1:length(data);
        if ~isnan(data(i));
            shotindex = i - m*floor(i/m);
            plot(s.xEvent(shotindex),s.yEvent(shotindex),'o','MarkerFaceColor',col(i,:),'MarkerEdgeColor',col(i,:),'Markersize',8)
        end
    end
    plot(s.xStation(ista),s.yStation(ista),'^k','MarkerFaceColor','k')
    axis('equal');
    ylabel('Y, km');
    xlabel('X, km');
    title(['Stations ' cell2mat(s.srStation.name(ista(1))) ' through ' cell2mat(s.srStation.name(ista(length(ista))))],'FontSize',11);
    caxis(lim);
    hcb = colorbar;
    h = gca;
    set(gcf,'currentaxes',hcb);
    if doNorm
        title('Norm Resid');
    else
        title('Residual, s');
    end
    set(gcf,'currentaxes',h);
end
