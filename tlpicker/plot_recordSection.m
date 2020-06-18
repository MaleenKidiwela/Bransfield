function [status,hRS,x0Trace,pickWidth,absScale] = plot_recordSection(tD, tMD, mP)
%Function to plot a record section for tlPicker
%
%Usage
%  [status,hRS,x0Trace,pickWidth,scale] = plot_recordSection(tD, tMD, mP)
%
%Inputs
%  tD  - Column matrix of traces
%  tMD - Structure of trace meta data with the following fields
%        Fields are described in get_segy.m
%  mP  - Structure of menu parameters with mP.plot a structure with the
%        following fields
%    xOption         - String indicated parameter to plot on X-axis
%                      It can be any of the vector fields of the tMD 
%                      (see get_segy.m) or 'sort' to indicate that the traces
%                      are plotted against their sort_index
%                      Note that the station option is not yet implemented
%                     
%    xlim0           - Lower x value of plot
%                      NaN (or equal to xlim1 indicates that it is set by data
%    xlim1           - Upper x value of plot
%                      NaN (or equal to xlim0 indicates that it is set by data
%    xScaleOption    - Single character to indicate type of trace scaling
%      'f' - Uniform scaling (also 'F')
%      'x' - Range dependent scaling (also 'X')
%      'u' - Uniform scaling to achieve desired maximum amplitude
%      'U' - Uniform scaling to achieve desired maximum amplitude
%                  including data outside plot area
%      'e' - Equal maximum amplitude
%      'E' - Equal maximum amplitude including data outside plot area
%      'r' - Range dependent scaling to achieve desired maximum amplitude
%      'R' - Range dependent scaling scaling to achieve desired maximum amplitude
%               including data outside plot area
%     xScale         - Trace scaling value
%     clip           - Clip traces with a deflection exceeding this value
%     wiggleOption   - Controls wiggle filling
%                      0 - No wiggle filling, line only
%                      1 - Wiggle filling on right hand side with line
%                      2 - Wiggle filling on right hand side without line
%                      3 - Blue-Red Color wiggle filling 
%                      4 - Red-Blue Color wiggle filling 
%                      -1, -2 - Wiggles filled on left hand side
%                      -3, -4 - As for 3,4 but traces plotted in opposite
%                      order
%     tlim0          - Lower time value of plotClip 
%                      NaN (or equal to tlim1 indicates that it is set by data)
%     tlim1          - Upper time value of plot
%                      NaN (or equal to tlim0 indicates that it is set by data)
%     redVel         - Reduction velocity for plotting (km/s)
%     demedian       - Prior to plotting the full loaded trace is demeaned
%                      This logical indicates whether the plotted portion of
%                      the trace is to have a median value of 0
%     defaultTitle   - Logical to indicate a default title is used which is
%                      given by
%                      ['{date} Station {station}; Channel {channel}; ' ...
%                       'Filter {filter}; Scale {scale}; Clip {clip}'];
%     title          - Title string which is interpreted as a text string
%                      except that
%                      {date} is replaced by the date in '26-Jan-2011' format
%                      {station} is replaced by the menu Station(s)
%                      {channel} is replaced by the menu Channel(s)
%                      {eventid} is replaced by the menu Events
%                      {filter} is replaced by the filter parameters in the format
%                      'Minimum phase, 4th order, [5 40] Hz'
%                      {scale} is replaced by the scale factor and Option
%                      in the fomat '0.0001f'
%                      {clip} is replaced by the clip parameter
%                      {xoption} is replaced by the xOption parameter and
%                      the sort parameters if XOption is 'sort' in format
%                      'sort w/ [range; channel]'
%                      @ is replaced by a carriage return
%     labelIncrement - Increment for trace labels
%     labelOption    - Up to three trace labels which are any of the vector
%                      fields of tMD (traceMetaData) - see get_segy.m
%     A few fields of mP.segy, mP.filter & mp.sort are used for the
%     {} options in title construction
%
%Outputs
%    status    - Status of execution
%                0 - okay
%                1 - Could not determine plotted center x value for traces
%    hRS       - Structure of handles (not implented)
%    x0Trace   - Vector of X values of traces corresponding to the input
%                traces.  NaN's for those traces outside plot area
%    pickWidth - Median spacing of traces (used to set plotted pick width)
%    absScale  - Vector with absolute scaling of each trace (0 if outside
%                plot limits)

%% Hardwired wiggle colors
wiggleColor1 = [1 0.2 0.2];
wiggleColor2 = [0.2 0.2 1];

%% Default output values
status = 0;
hRS = [];
x0Trace = [];
pickWidth = [];
absScale = [];

%% Process Inputs
if ~isfield(tMD,'static')
  tMD.static = zeros(tMD.ntrace,1);
end
if ~isfield(tMD,'mute')  
  tMD.mute = Inf(tMD.ntrace,1);
end  
  
%% Get center X values for traces
if strcmpi(deblank(mP.plot.xOption), 'sort')
  x0(tMD.sortIndex) = 1:tMD.ntrace;
else
  names = fieldnames(tMD);
  k = find(strcmpi(names, mP.plot.xOption));
  if isempty(k)
    status = 1;
    return
  else
    eval(['x0 = double(tMD.' names{k} '(:))'';']);
    if length(x0) ~= tMD.ntrace
      status = 1;
      return
    end
  end
end

% Create outputs before reordering
x0Trace = x0;
pickWidth = median(diff(sort(x0Trace)));
absScale = zeros(tMD.ntrace,1);
  
%% Find traces with X values within X limits
xlimits = zeros(1,2);
if ~isnan(mP.plot.xlim0)
  xlimits(1) = mP.plot.xlim0;
else
  xlimits(1) = -Inf;
end
if ~isnan(mP.plot.xlim1)
  xlimits(2) = mP.plot.xlim1;
else
  xlimits(2) = Inf;
end
if ~diff(xlimits)
  xlimits = [ -Inf Inf];
end
if diff(xlimits)>0
  indexPlot = x0>=xlimits(1) & x0<=xlimits(2);
else
  indexPlot = x0>=xlimits(2) & x0<=xlimits(1);
end  
x0Trace(~indexPlot) = NaN;
indexPlot = find(indexPlot);
x0 = x0(indexPlot);
tD = tD(:,indexPlot);
tMD = subset_oneElementStructure(tMD,'ntrace',indexPlot);
if ~tMD.ntrace
  status = 2;
  return
end

%% Sort Traces
[x0,i] = sort(x0);
reverseX = rem((mP.plot.wiggleOption>0) + (diff(xlimits)<0),2);
if reverseX
  i = fliplr(i);
  x0 = fliplr(x0);
end
indexPlot = indexPlot(i);
tD = tD(:,i);
tMD = subset_oneElementStructure(tMD,'ntrace',i);

%% Determine time limits
tlimits = zeros(1,2);
if ~isnan(mP.plot.tlim0)
  tlimits(1) = mP.plot.tlim0;
else
  tlimits(1) = -Inf;
end
if ~isnan(mP.plot.tlim1)
  tlimits(2) = mP.plot.tlim1;
else
  tlimits(2) = Inf;
end
if ~diff(tlimits)
  tlimits = [ -Inf Inf];
end
reverseT = tlimits(1) > tlimits(2);

%% Create matrices of sample values and times
% maxAmpPlot and maxAmpAll store trace maximum absolute amplitudes for 
maxAmpPlot = zeros(1,tMD.ntrace);
maxAmpAll = zeros(1,tMD.ntrace);
nsampPlot = zeros(1,tMD.ntrace);
XMAT = [];
TMAT = [];
for i = 1:tMD.ntrace
  % Time
  tvector = tMD.t0(i) - tMD.static(i) + (0:tMD.nsamp(i)-1)'/tMD.samprate(i);
  tmute = tMD.mute(i) - tMD.static(i);
  if ~isnan(mP.plot.redVel)
    if mP.plot.redVel
      tvector = tvector - abs(tMD.range(i))./(mP.plot.redVel);
      tmute = tmute - abs(tMD.range(i))./(mP.plot.redVel);
    end
  end
  if ~reverseT
    iplot = find(tvector>=tlimits(1) & tvector<=tlimits(2));
  else
    iplot = find(tvector>=tlimits(2) & tvector<=tlimits(1));
  end
  tvector = tvector(iplot);
  ntvector = length(tvector);
  % Values
  xvector = tD(1:tMD.nsamp(i),i);
  xvector = xvector - mean(xvector(~isnan(xvector)));
  maxAmpAll(i) = max(xvector(~isnan(xvector)));
  xvector = xvector(iplot);
  % Color options require zero crossings to be added
  if abs(mP.plot.wiggleOption) > 2 
    i0 = find(diff(sign(xvector)));
    ic = i0 - xvector(i0) ./ (xvector(i0+1) - xvector(i0));
    tcross = tvector(1) + (ic-1)/tMD.samprate(i);
    tvector = [tvector; tcross; tvector(1)-1e-10; tvector(end)+1e-10];
    xvector = [xvector; zeros(length(ic)+2,1)];
    [tvector,is] = sort(tvector);
    xvector = xvector(is);
    ntvector = length(tvector);
  end
  % Muting  
  j = find(tvector<tmute,1,'last');
  if ~isempty(j);
    xvector(j+1:end) = NaN;
    if mP.plot.demedian
      xvector = xvector-median(xvector(~isnan(xvector)));
    end
    if any(~isnan(xvector))
      maxAmpPlot(i) = max(abs(xvector(~isnan(xvector))));
    else
      maxAmpPlot(i) = 1;
    end
    nsampPlot(i) = j; 
  else
    xvector(1:end) = NaN;
    maxAmpPlot(i) = 0;
    nsampPlot(i) = 0;    
  end
  if i==1
    XMAT = NaN(ntvector,tMD.ntrace);
    TMAT = NaN(ntvector,tMD.ntrace);
  end
  TMAT(1:ntvector,i) = tvector(:);
  XMAT(1:ntvector,i) = xvector;
end
i = max(nsampPlot);
XMAT = XMAT(1:i,:);
TMAT = TMAT(1:i,:);
if isempty(XMAT)
  status = 3;
  return
end
  
%% Scaling with the following options
absScale = zeros(tMD.ntrace,1);
scale = mP.plot.xScale;
if strcmp(mP.plot.xScaleOption,'u')
  scale = scale / max(maxAmpPlot);
elseif strcmp(mP.plot.xScaleOption,'U')
  scale = scale / max(maxAmpAll);
elseif strcmp(mP.plot.xScaleOption,'r')
  scale = scale / max(maxAmpPlot.*abs(tMD.range'));
elseif strcmp(mP.plot.xScaleOption,'R')
  scale = scale / max(maxAmpAll.*abs(tMD.range'));
end
for i = 1:tMD.ntrace
  if strcmpi(mP.plot.xScaleOption,'f')
    absScale(indexPlot(i)) = scale;
  elseif strcmpi(mP.plot.xScaleOption,'x')
    absScale(indexPlot(i)) = scale*abs(tMD.range(i));    
  elseif strcmp(mP.plot.xScaleOption,'u')
    absScale(indexPlot(i)) = scale;
  elseif strcmp(mP.plot.xScaleOption,'U')
    absScale(indexPlot(i)) = scale;
  elseif strcmp(mP.plot.xScaleOption,'r')
    absScale(indexPlot(i)) = scale*abs(tMD.range(i));
  elseif strcmp(mP.plot.xScaleOption,'R')
    absScale(indexPlot(i)) = scale*abs(tMD.range(i));
  elseif strcmp(mP.plot.xScaleOption,'e')
    absScale(indexPlot(i)) = scale/maxAmpPlot(i);
  elseif strcmp(mP.plot.xScaleOption,'E')
    absScale(indexPlot(i)) = scale/maxAmpAll(i);
  end
  XMAT(:,i) = XMAT(:,i)*absScale(indexPlot(i));
end

%% Clipping
if ~isnan(mP.plot.clip)
  if mP.plot.clip
    XMAT(XMAT>abs(mP.plot.clip)) = abs(mP.plot.clip);
    XMAT(XMAT<-abs(mP.plot.clip)) = -abs(mP.plot.clip);
  end
end

%% Add center X values to traces
XMAT = bsxfun(@plus, XMAT, x0);

%% Make Axes
if ~isinf(diff(xlimits))
  xlimits2 = sort(xlimits);
else
  xlimits2 = [min(min(XMAT)) max(max(XMAT))];
end
if ~isinf(diff(tlimits))
  tlimits2 = sort(tlimits);
else
  tlimits2 = [min(min(TMAT)) max(max(TMAT))];
end
clf
h1 = axes('xlim',xlimits2,'ylim',tlimits2);

if diff(xlimits)<0
  set(h1,'xdir','rev')
end
if diff(tlimits)<0
  set(h1,'ydir','rev')
end
if ~isempty(mP.plot.labelOption)
  pos=[.10 .08 .85 0.7750];
else
  pos=[.10 .08 .85 0.8150];
end
set(h1,'units','normalized','position',pos)


%% Plot the traces
if reverseX
  extrema = min(xlimits2);
else
  extrema = max(xlimits2);
end  
if ~isnan(mP.plot.wiggleOption)
  if mP.plot.wiggleOption
    for i=1:tMD.ntrace
      if abs(mP.plot.wiggleOption)<=2
        % B&W filled wiggles
        if nsampPlot(i)
          xp = [extrema; x0(i); x0(i);  XMAT(1:nsampPlot(i),i);  ...
              x0(i);  x0(i);  extrema];
          tp = [min(tlimits2);  min(tlimits2);  TMAT(1,i);  TMAT(1:nsampPlot(i),i);  ...
              TMAT(nsampPlot(i),i);  max(tlimits2);  max(tlimits2)];
          xp(isnan(xp)) = x0(i);
          patch(xp,tp,'black','edgecolor','none');
          patch([extrema;  x0(i);  x0(i);  extrema], ...
              [min(tlimits2);  min(tlimits2);  max(tlimits2);  max(tlimits2)], ...
              'white','edgecolor','none');
        end
      else
        % Color filled wiggles
        if abs(mP.plot.wiggleOption)==3
          patch(max(XMAT(1:nsampPlot(i),i),x0(i)),TMAT(1:nsampPlot(i),i),wiggleColor1,'edgecolor','none')
          patch(min(XMAT(1:nsampPlot(i),i),x0(i)),TMAT(1:nsampPlot(i),i),wiggleColor2,'edgecolor','none')
        elseif abs(mP.plot.wiggleOption)==4
          patch(max(XMAT(1:nsampPlot(i),i),x0(i)),TMAT(1:nsampPlot(i),i),wiggleColor2,'edgecolor','none')
          patch(min(XMAT(1:nsampPlot(i),i),x0(i)),TMAT(1:nsampPlot(i),i),wiggleColor1,'edgecolor','none') 
        end
      end
    end
  end
  if abs(mP.plot.wiggleOption)<2
    doline = true;
  else
    doline = false;
  end
else
  doline = true;
end
if doline
  for i=1:tMD.ntrace
    if nsampPlot(i)
      line(XMAT(1:nsampPlot(i),i),TMAT(1:nsampPlot(i),i),'color','black');
    end
  end
end

%% Create Labelled Axes (Station option needs work)
xtick = get(h1,'xtick');
set(h1,'xticklabel',num2str(xtick'));
if strcmpi(mP.plot.xOption,'sort')
  xlabel('Sorted Trace Index')
elseif strcmpi(mP.plot.xOption,'sortindex')
  xlabel('Sorted Trace Index')
elseif strcmpi(mP.plot.xOption,'range')
  xlabel('Range, km')
elseif ~isempty(findstr(lower(mP.plot.xOption),'lat'))
  xlabel('Latitude, \circ')
elseif ~isempty(findstr(lower(mP.plot.xOption),'lon'))
  xlabel('Longitude, \circ')
elseif findstr(lower(mP.plot.xOption),'x')
  if any(findstr(lower(mP.plot.xOption),'x')==1)
    xlabel('X, km')
  end
elseif findstr(lower(mP.plot.xOption),'y')
  if any(findstr(lower(mP.plot.xOption),'y')==1)
    xlabel('Y, km')
  end
elseif strcmpi(mP.plot.xOption,'eventid')
  xlabel('Event ID')
elseif strcmpi(mP.plot.xOption,'eventindex')
  xlabel('Event Index')
elseif strcmpi(mP.plot.xOption,'stationindex')
  xlabel('Station Index')
elseif strcmpi(mP.plot.xOption,'channel')
  xlabel('Channel')
elseif strcmpi(mP.plot.xOption,'azimuth')
  xlabel('Azimuth, \circ')
elseif strcmpi(mP.plot.xOption,'station')
  xlabel('Station')
end
if ~isnan(mP.plot.redVel)
  if mP.plot.redVel
     ylabel(['Time - Range/' num2str(mP.plot.redVel) ', s'], ...
                    'fontsize',get(gcf,'defaultaxesfontsize')+2);
  else
     ylabel('Time, s','fontsize',get(gcf,'defaultaxesfontsize')+2);
  end
else
  ylabel('Time, s','fontsize',get(gcf,'defaultaxesfontsize')+2);
end

%% Plot Trace labels
if ~isempty(mP.plot.labelOption) && ~isempty(mP.plot.labelIncrement)
  if mP.plot.labelIncrement
    if reverseT
      tLabel = tlimits2(1) - 0.01*diff(tlimits2);
    else
      tLabel = tlimits2(2) + 0.01*diff(tlimits2);
    end  
    labels = lower(string2cell(mP.plot.labelOption));
    labels = strrep(labels,'long','lon');
    labels = strrep(labels,'index','Index');
    j = 0;
    values = nan(tMD.ntrace,length(labels));
    for i = 1:length(labels);
      j = j + 1;
      try
        eval(['values(:,j) = tMD.' labels{i} ';'])
        if ~isempty(strfind(labels{i},'x')) || ...
           ~isempty(strfind(labels{i},'y')) || ...
           ~isempty(strfind(labels{i},'range'))
          units(j) = {' km'};
        elseif strcmp(labels{i},'azimuth')
          units(j) = {'\circ'};
        else
          units(j) = {''};
        end
      catch
        warning('Invalid label option')
        j = j - 1;
      end
    end
    nvalues = j;
    values = round(values*10)/10;
    clear string
    for i = mP.plot.labelIncrement:mP.plot.labelIncrement:tMD.ntrace
      for j = 1:nvalues  
        string(j) = cellstr([num2str(values(i,j)) units{j}]);
      end
      text(x0(i),tLabel,string, ...
          'verticalalignment','middle','rotation',90, ...
          'horizontalalignment','left','fontsize',10);

    end
  end
end

%% Plot title
if mP.plot.defaultTitle
  string = '{date} Station {station}; Channel {channel}; Filter {filter}; Scale {scale}; Clip {clip}';
else
  string = mP.plot.title;
end
% Date 
index = fliplr(strfind(lower(string),'{date}'));
for i = index;
  string = [string(1:i-1) date string(i+6:end)];
end
% Station
index = fliplr(strfind(lower(string),'{station}'));
for i = index;
  string = [string(1:i-1) mP.segy.station string(i+9:end)];
end
% Channel
index = fliplr(strfind(lower(string),'{channel}'));
for i = index;
  string = [string(1:i-1) num2str(mP.segy.channel) string(i+9:end)];
end
% Filter
index = fliplr(strfind(lower(string),'{filter}'));
for i = index;
  string1 = '';
  if mP.filter.zeroPhase
    string1 = 'Zero phase, ';
  else
    string1 = 'Minimum phase, ';
  end
  if isempty(mP.filter.order) || ~mP.filter.order
    string2 = '';
  elseif mP.filter.order == 1
    string2 = '1st';
  elseif mP.filter.order == 2
    string2 = '2nd';
  elseif mP.filter.order == 3
    string2 = '3rd';
  elseif mP.filter.order == 4
    string2 = [num2str(mP.filter.order) 'th'];
  end
  string = [string(1:i-1) string2 ' order, ' string1 ' [' ...
           num2str([mP.filter.lim0 mP.filter.lim1]) ...
           '] Hz' string(i+8:end)];
end
% Event ID
index = fliplr(strfind(lower(string),'{eventid}'));
for i = index;
  string = [string(1:i-1) num2str(mP.segy.eventid) string(i+9:end)];
end   
% Scale
index = fliplr(strfind(lower(string),'{scale}'));
for i = index;
  string = [string(1:i-1) num2str(mP.plot.xScale) ...
            mP.plot.xScaleOption string(i+7:end)];
end   
index = fliplr(strfind(lower(string),'{clip}'));
for i = index;
  string = [string(1:i-1) num2str(mP.plot.clip) string(i+6:end)];
end   
index = fliplr(strfind(lower(string),'{xoption}'));
if ~isempty(index)
  if strcmpi(mP.plot.xOption,'sort')
    string1 = 'sort w/ [';
    if ~isempty(mP.sort(1).name)
      string1 = [string1 mP.sort(1).name];
    end
    if ~isempty(mP.sort(2).name)
      string1 = [string1 '; ' mP.sort(2).name];
    end
    if ~isempty(mP.sort(3).name)
      string1 = [string1 '; ' mP.sort(3).name];
    end
    if ~isempty(mP.sort(4).name)
      string1 = [string1 '; ' mP.sort(4).name];
    end
    string1 = [string1 ']'];
  else
    string1 = mP.plot.xOption;
  end
  for i = index;
    string = [string(1:i-1) num2str(mP.plot.xOption) string(i+6:end)];
  end
end
string = string2cell(string,'@');
if isempty(mP.plot.labelOption)
  text(mean(xlim),max(ylim)+abs(diff(ylim))*.05,string,'horizontalalignment','center','fontsize', ...
         get(gcf,'defaultaxesfontsize')+2,'interpreter','none');
else
  text(mean(xlim),max(ylim)+abs(diff(ylim))*.15,string,'horizontalalignment','center','fontsize', ...
         get(gcf,'defaultaxesfontsize')+2,'interpreter','none');
end
