function [traceMetaData,traceData] = ...
  get_segy(paramSegy,segyCatalog,srEvent,srStation,srGeometry,srStationRotation)
% Reads cataloged segy data for stingray/tomolab applications
%
% Usage
%   [traceMetaData,traceData] = ...
%   get_segy(paramSegy,segyCatalog,srEvent,srStation,srGeometry,srStationRotation)
% 
% Inputs
%   paramSegy - Structure determining which traces to load with fields
%     station     - Station name(s) as a string 
%     eventid     - Event ID(s) as a string
%                     Special ETOMO
%                     -1 to -45 = lines 1 to 45
%                     -52 stacked eastern undershoot line (stack created on the fly)
%                     -55 stacked western undershoot line (stack created on the fly) 
%                     -53 stacked line inboard of easter undershoot (created on fly)
%                     -62 stacked eastern undershoot line (stack exists as SEGY)
%                     -65 stacked western undershoot line (stack exists as SEGY)
%                     -63 stacked line in board of eastern undershoot (stack exists as SEGY)
%     channel     - Channel number(s) as a string
%                     1, 2, 3, 4 = Vert, Horiz1, Horiz2, Hydrophone
%                     5, #   = Vertical + # * Hydrophone
%                     6, 7   = Radial and Transverse channels
%                     -6, -7 = Radial and Transverse channels with
%                               alternate rotation (based on differential
%                               scaling of channels 2 and 3
%     excludeBad  - Logical to ignore traces tagged as bad in catalog
%     tlim0       - Starting time of data to load
%     tlim1       - Ending time of data to load
%     redVel      - Reduction velocity (km/s)
%     select      - Structure vector with additional selection criteria.  
%                   Each element of the structure has the following fields
%       name - Name of the selection parameter.  Current options are
%              'range'             - Range (km)
%              'azimuth'           - Azimuth (degrees)
%              'xmid'              - X value of source-receiver mid point
%              'ymid'              - Y value ...
%              'latmid'            - Latitude ...
%              'lonmid', 'longmid' - Longitue ...
%              'xsta'              - X value of receiver location
%              'ysta'              - Y value ...
%              'latsta'            - Latitude ...
%              'lonsta', 'longsta' - Longitue ...
%              'xevt'              - X value of source location
%              'yevt'              - Y value ...
%              'latevt'            - Latitude ...
%              'lonevt', 'longsta' - Longitue ...          
%       lim0 - Lower limit of selection parameter
%       lim1 - upper limit of selection paramter
%     pickSelect  - Stucture to select only traces with picks
%                   Picks maybe input in 3 formats
%                   1. tlArrival structure/file and phase(s)
%                   2. tomoLab predictions via a tomoLab inversion
%                      directory, phase(s) and iteration
%                   3. tlPicker pick flat ascii pick files via the
%                      directory and phases(s)
%                   The structure has 3 fields
%       file      - File (or directory or structure) name
%       phase     - String of phase(s)
%       iteration - Iteration number (option 2 only)
%  segyCatalog - structure created by creat_segyCatalog_@@@@.m
%                (where @@@ is the experiment) that indexes the SEGY data
%                for fast searching and reading.  
%                It contains the following fields:
%     Vectors dimensioned by the number of SEGY files:
%       filename     - Cell vector of N SEGY filenames
%       format       - Binary header data sample format code for each of 
%                      N SEGY files
%                      (see: http://www.seg.org/SEGportalWEBproject/prod/
%                      SEG-Publications/Pub-Technical-Standards/Documents/
%                      seg_y_rev1.pdf)
%       littleEndian - Logical to indicate little or big Endian file
%                      formats for each of N SEGY files
%     3-D matrix dimensioned by nchan by srEvent.nevt by srStation.nsta*
%       fileIndex    - Index of file containing each trace
%       byte0        - Byte index of 1st sample of trace in file
%       nsamp        - Number of samples in trace
%       samprate     - Sample rate of trace (Hz)
%       delay        - Delay time for 1st sample (s)
%       bad          - Trace is flagged as bad
%       xmid         - X mid point (km) of trace
%       ymid         - Y mid point (km) of trace
%       range        - Source-receiver range of trace (km)
%       azimuth      - Source-receiver azimuth of trace (degrees)
%       * All matrices are collapsed along invarient dimensions 
%  srEvent     - stingray event structure
%  srStation   - stingray station structure
%  srGeometry  - stingray geometry structure
%     
% Outputs
%   traceData     - Matrix of seismic traces in columns (NaNs for no data)
%   traceMetaData - Scalar stucture of vectors of trace meta data
%     ntrace        - Number of traces
%     nsamp         - Number of samples for each trace
%     samprate      - Sample rate (Hz)
%     t0            - Time (s) of 1st sample (relative to event time)
%     eventid       - Event ID
%     eventIndex    - Event index within srEvent
%     station       - Station name
%     stationList   - List of stations with at least one trace
%     stationIndex  - Station index within srStation
%     channel       - Channel number
%     range         - Range (km)
%     azimuth       - Azimuth (degrees)
%     xmid          - Source-receiver mid point X value (km)
%     ymid          - Source-receiver mid point Y value (km)
%     xsta          - Receiver X value (km)
%     ysta          - Receiver Y value (km)
%     xevt          - Source X value (km)
%     yevt          - Source Y value (km)
%     lonmid        - Source-receiver mid point longitude (degrees)
%     latmid        - Source-receiver mid point latitude (degrees)
%     lonsta        - Source longitude (degrees)    
%     latsta        - Source latitude (degrees)
%     lonevt        - Receiver longitude (degrees)
%     latevt        - Receiver latitude (degrees)
%     nstack        - Number of traces stacked internally 
%
%     Three fields use for plotting are initially set to null values
%     sortIndex     - Sorting index of traces is set to 1:ntraces
%     static        - Static correction vector (s) is set to 0
%     mute          - Mute time vector (s) is set to Inf
%     nstack        - Number of traces stacked internally
%
% 4/15/11 - Fixed bug in terms of accessing data because single precision
% can be insufficient for fseek and it got set to single precision because
% i0 was singe precision.  Also ensure all numerical fields in
% traceMetaData are double precision

MAXTRACE = 6000;

[nchan,nevt,nsta] = size(segyCatalog.byte0);

% Selection starts with all available traces
select = segyCatalog.byte0>0;

% Select stations
discardNow = true(nsta,1);
name = string2cell(paramSegy.station);
if ~any(strcmp(name,'0'))
  for i = 1:length(name)
    discardNow(strcmp(srStation.name,name(i))) = false;
  end
  if all(discardNow)
    warning('Failure in get_segy.m - No valid stations specified')
    traceMetaData = [];
    traceData = [];
    return  
  else
    select(:,:,discardNow) = false;
  end
end

% Select eventids
discardNow = true(nevt,1);
vector = round(string2vector(paramSegy.eventid));
vector2 = [];
% ETOMO Specific Start
warnMissedEvent = true;
for i = 1:length(vector)
  if (vector(i)<0 && vector(i)>=-45) || vector(i)==-62 || vector(i)==-65 || vector(i)==-63
    warnMissedEvent = false;
    vector2 = [vector2 -vector(i)*1000+(0:999)];
  
  elseif vector(i) == -52
    warnMissedEvent = false;
    vector2 = [vector2 7000+(0:999) 8000+(0:999) 39000+(0:999)];
  
  elseif vector(i) == -53
    warnMissedEvent = false;
    vector2 = [vector2 4000+(0:999) 6000+(0:999)];
  
  elseif vector(i) == -55
    warnMissedEvent = false;
    vector2 = [vector2 1000+(0:999) 42000+(0:999)];
  else
    vector2 = [vector2 vector(i)];
  end
end
vector = vector2;
% ETOMO Specific End

if ~any(vector==0)
  tf = vector_match(srEvent.id,vector);
  if sum(tf)~=length(vector)
    if warnMissedEvent
      warning('Some values of eventid not found in get_segy.m')
    end
  end
  discardNow(tf) = false;
  if all(discardNow)
    warning('Failure in get_segy.m - No valid events specified')
    traceMetaData = [];
    traceData = [];
    return  
  else
    select(:,discardNow,:) = false;
  end
end

% Select channels 
% ETOMO Specific Start
% (5 is a weighted sum of channels 1 and 4 with the weight of channel 4 as the next number in string)
% (6 = radial & 7 = tangential with uniform scaling of rotated channels 2/3)
% (-6 and -7 = as for above except with relative scaling or rotated channels)
channel5 = false;
vector = string2vector(paramSegy.channel);
if ~any(vector==0)
  discardNow = true(nchan,1);
  discardAtEnd = true(nchan,1);
  for i = 1:length(vector)
    if isnan(channel5)
      channel5 = vector(i);
    elseif vector(i)>=1 && vector(i)<=4
      discardNow(vector(i)) = false;
      discardAtEnd(vector(i)) = false;
    elseif vector(i)==5 
      discardNow([1 4]) = false;
      if i == length(vector)
        warning('No channel 4 scaling value provided following channel 5 selection')
        discardNow = discardNow | true;
        break
      else
        channel5 = NaN;
      end
    elseif vector(i)==6 || vector(i)==7 || vector(i)==-6 || vector(i)==-7
      discardNow([2 3]) = false;
    end
  end
  if all(discardNow)
    warning('Failure in get_segy.m - No valid channels specified')
    traceMetaData = [];
    traceData = [];
    return  
  else
    select(discardNow,:,:) = false;
  end
end
% ETOMO Specific End

% Reject bad traces
if paramSegy.excludeBad
  select = select & ~segyCatalog.bad;
end

% Convert mid-point x,y to longitude and latitude if necessary
if srGeometry.tf_latlon
  for isel = 1:length(paramSegy.select)
    if strcmpi(paramSegy.select(isel).name,'lonmid') || ...
       strcmpi(paramSegy.select(isel).name,'latmid')
      [lonmid,latmid] = xy2map(segyCatalog.xmid,segyCatalog.ymid,srGeometry);
      break
    end
  end
end

% Select based on selection parameters
for isel = 1:length(paramSegy.select)
  if ~isempty(paramSegy.select(isel).name) && ...
    (~isempty(paramSegy.select(isel).lim0) || ~isempty(paramSegy.select(isel).lim1))
    switch lower(paramSegy.select(isel).name)

      case 'range'
        if ~isempty(paramSegy.select(isel).lim0)
          select = select & segyCatalog.range>=paramSegy.select(isel).lim0;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          select = select & segyCatalog.range<=paramSegy.select(isel).lim1;
        end

      case 'azimuth'
        if ~isempty(paramSegy.select(isel).lim0) && ~isempty(paramSegy.select(isel).lim1)
          if paramSegy.select(isel).lim0 <= paramSegy.select(isel).lim1
            select = select & segyCatalog.azimuth>=paramSegy.select(isel).lim0 ...
                            & segyCatalog.azimuth<=paramSegy.select(isel).lim1;
          else
            select = select & (  segyCatalog.azimuth>=paramSegy.select(isel).lim0 ...
                               | segyCatalog.azimuth<=paramSegy.select(isel).lim1);
          end
        else
          if ~isempty(paramSegy.select(isel).lim0)
            select = select & segyCatalog.azimuth>=paramSegy.select(isel).lim0;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            select = select & segyCatalog.azimuth<=paramSegy.select(isel).lim1;
          end  
        end

      case 'xmid'
        if ~isempty(paramSegy.select(isel).lim0)
          select = select & segyCatalog.xmid>=paramSegy.select(isel).lim0;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          select = select & segyCatalog.xmid<=paramSegy.select(isel).lim1;
        end

      case 'ymid'
        if ~isempty(paramSegy.select(isel).lim0)
          select = select & segyCatalog.ymid>=paramSegy.select(isel).lim0;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          select = select & segyCatalog.ymid<=paramSegy.select(isel).lim1;
        end

      case {'lonmid', 'longmid'}
        if ~srGeometry.tf_latlon
          warning('Longitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            select = select & lonmid>=paramSegy.select(isel).lim0;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            select = select & lonmid<=paramSegy.select(isel).lim1;
          end
        end
        
      case 'latmid'
        if ~srGeometry.tf_latlon
          warning('Latitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            select = select & latmid>=paramSegy.select(isel).lim0;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            select = select & latmid<=paramSegy.select(isel).lim1;
          end
        end
        
      case 'xsta'
        if ~isempty(paramSegy.select(isel).lim0)
          discardNow = srStation.x<paramSegy.select(isel).lim0;
          select(:,:,discardNow) = false;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          discardNow = srStation.x>paramSegy.select(isel).lim1;
          select(:,:,discardNow) = false;
        end
        
      case 'ysta'
        if ~isempty(paramSegy.select(isel).lim0)
          discardNow = srStation.y<paramSegy.select(isel).lim0;
          select(:,:,discardNow) = false;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          discardNow = srStation.y>paramSegy.select(isel).lim1;
          select(:,:,discardNow) = false;
        end
        
      case {'lonsta', 'longsta'}
        if srGeometry.tf_latlon==false
          warning('Longitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            discardNow = srStation.longitude<paramSegy.select(isel).lim0;
            select(:,:,discardNow) = false;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            discardNow = srStation.longitude>paramSegy.select(isel).lim1;
            select(:,:,discardNow) = false;
          end
        end
        
      case 'latsta'
        if srGeometry.tf_latlon==false
          warning('Latitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            discardNow = srStation.latitude<paramSegy.select(isel).lim0;
            select(:,:,discardNow) = false;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            discardNow = srStation.latitude>paramSegy.select(isel).lim1;
            select(:,:,discardNow) = false;
          end
        end
        
      case 'xevt'
        if ~isempty(paramSegy.select(isel).lim0)
          discardNow = srEvent.x<paramSegy.select(isel).lim0;
          select(:,discardNow,:) = false;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          discardNow = srEvent.x>paramSegy.select(isel).lim1;
          select(:,discardNow,:) = false;
        end
        
      case 'yevt'
        if ~isempty(paramSegy.select(isel).lim0)
          discardNow = srEvent.y<paramSegy.select(isel).lim0;
          select(:,discardNow,:) = false;
        end
        if ~isempty(paramSegy.select(isel).lim1)
          discardNow = srEvent.y>paramSegy.select(isel).lim1;
          select(:,discardNow,:) = false;
        end
                
      case {'lonevt', 'longevt'}
        if srGeometry.tf_latlon==false
          warning('Longitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            discardNow = srEvent.longitude<paramSegy.select(isel).lim0;
            select(:,discardNow,:) = false;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            discardNow = srEvent.longitude>paramSegy.select(isel).lim1;
            select(:,discardNow,:) = false;
          end
        end
        
      case 'latevt'
        if srGeometry.tf_latlon==false
          warning('Latitude limits ignored in get_segy.m for easting/northing geometry')
        else
          if ~isempty(paramSegy.select(isel).lim0)
            discardNow = srEvent.latitude<paramSegy.select(isel).lim0;
            select(:,discardNow,:) = false;
          end
          if ~isempty(paramSegy.select(isel).lim1)
            discardNow = srEvent.latitude>paramSegy.select(isel).lim1;
            select(:,discardNow,:) = false;
          end
        end
                
      otherwise
        warning(['Unrecognized Parameter name ' ...
               paramSegy.select(isel).name ' in get_segy.m'])
        
    end
  end
end

% Indicies of presently selected traces
icat = find(select);
[channel,ievt,ista] = ind2sub([nchan nevt nsta],icat);

% Apply a pick file to the selected traces
if ~isempty(paramSegy.pickSelect)
  if ~isempty(paramSegy.pickSelect.DFS)
    if paramSegy.pickSelect.DFS(end) == '/'
      paramSegy.pickSelect.DFS = paramSegy.pickSelect.DFS(1:end-1);
    end
    if ~exist(paramSegy.pickSelect.DFS)
      warning('Pick directory/file/structure for trace selection does not exist');
    else
      tMDtemp.station = srStation.name(ista);
      tMDtemp.eventid = srEvent.id(ievt);
      tMDtemp.channel = channel(:);
      [pick, status] = get_pickMatrix(paramSegy.pickSelect.DFS, ...
                       paramSegy.pickSelect.phase, ...
                       paramSegy.pickSelect.chanIter, tMDtemp);
      if status==2
        warning('get_segy: DFS is not a valid directory/file/structure - Ignored')
      elseif status==3
        warning('get_segy: DFS contains multiple picks for a station/eventid/channel/phase - Ignored')
      elseif status==4
        warning('get_segy: DFS needs a phase for a tlPick file - Ignored')
      else
        if size(pick,2) == 1
          select = ~isnan(pick);
        else
          select = any(~isnan(pick'))';
        end
        icat = icat(select);
        channel = channel(select);
        ievt = ievt(select);
        ista = ista(select);
      end
    end
  end
end

% Check number of traces
ntrace = length(icat);
if ~ntrace
  warning('get_segy.m - No traceMetaDatas satisfy selection criterion')
  traceMetaData = [];
  traceData = [];
  return  
elseif ntrace>MAXTRACE
  warning('get_segy.m - Too many traces - Increase MAXTRACE parameter in this function')
  traceMetaData = [];
  traceData = [];
  return
end

% Get catalog data for selected traceMetaDatas
fileIndex = zeros(ntrace,1);
byte0 = zeros(ntrace,1);
nsamp = zeros(ntrace,1);
samprate = zeros(ntrace,1);
delay = zeros(ntrace,1);
for i = 1:ntrace
  fileIndex(i) = segyCatalog.fileIndex( min(channel(i),size(segyCatalog.fileIndex,1)), ...
                                        min(ievt(i),size(segyCatalog.fileIndex,2)), ...
                                        min(ista(i),size(segyCatalog.fileIndex,3)));
  byte0(i) = segyCatalog.byte0( min(channel(i),size(segyCatalog.byte0,1)), ...
                                min(ievt(i),size(segyCatalog.byte0,2)), ...
                                min(ista(i),size(segyCatalog.byte0,3)));
  nsamp(i) = segyCatalog.nsamp( min(channel(i),size(segyCatalog.nsamp,1)), ...
                                min(ievt(i),size(segyCatalog.nsamp,2)), ...
                                min(ista(i),size(segyCatalog.nsamp,3)));
  samprate(i) = segyCatalog.samprate( min(channel(i),size(segyCatalog.samprate,1)), ...
                                      min(ievt(i),size(segyCatalog.samprate,2)), ...
                                      min(ista(i),size(segyCatalog.samprate,3)));
  delay(i) = segyCatalog.delay( min(channel(i),size(segyCatalog.delay,1)), ...
                                min(ievt(i),size(segyCatalog.delay,2)), ...
                                min(ista(i),size(segyCatalog.delay,3)));
end
range = double(segyCatalog.range(icat));
azimuth = double(segyCatalog.azimuth(icat));
xmid = double(segyCatalog.xmid(icat));
ymid = double(segyCatalog.ymid(icat));
  
% Sort by fileIndex for efficient reading 
[fileIndex,index] = sort(fileIndex);
byte0 = byte0(index);
nsamp = nsamp(index);
samprate = samprate(index);
delay = delay(index);
channel = channel(index);
ievt = ievt(index);
ista = ista(index);
range = range(index);
azimuth = azimuth(index);
xmid = xmid(index);
ymid = ymid(index);
% mute = mute(index);
% static = static(index);

% Reduction slowness
reductionSlowness = 0;
if ~isempty(paramSegy.redVel)
  if paramSegy.redVel
    reductionSlowness = 1/paramSegy.redVel;
  end
end

% Indicies of samples to read
if ~isempty(paramSegy.tlim0)
  i0 = ceil((paramSegy.tlim0 + range*reductionSlowness - delay) .*samprate)+1;
  i0 = max(1,min(i0,nsamp));
else
  i0 = ones(ntrace,1);
end
if ~isempty(paramSegy.tlim1)
  i1 = ceil((paramSegy.tlim1 + range*reductionSlowness - delay) .*samprate);
  i1 = max(1,min(i1,nsamp));
else
  i1 = nsamp;
end
i1 = max(i1,i0);
nsamp = (i1-i0)+1;
t0 = delay + (i0-1)./samprate;

% Assign traceMetaData metadata
traceData = NaN(max(nsamp),ntrace);
traceMetaData.ntrace = ntrace;
traceMetaData.nsamp = nsamp;
traceMetaData.samprate = samprate;
traceMetaData.t0 = t0;
traceMetaData.eventid = srEvent.id(ievt);
traceMetaData.eventIndex = ievt;
traceMetaData.station = srStation.name(ista);
traceMetaData.stationList = unique(traceMetaData.station);
traceMetaData.stationIndex = ista;
traceMetaData.channel = channel;
traceMetaData.range = range;
traceMetaData.azimuth = azimuth;
traceMetaData.xmid = xmid;
traceMetaData.ymid = ymid;
traceMetaData.xsta = srStation.x(ista);
traceMetaData.ysta = srStation.y(ista);
traceMetaData.xevt = srEvent.x(ievt);
traceMetaData.yevt = srEvent.y(ievt);
if srGeometry.tf_latlon==true
  [traceMetaData.lonmid,traceMetaData.latmid] = xy2map(traceMetaData.xmid,traceMetaData.ymid,srGeometry);
  traceMetaData.lonsta = srStation.longitude(ista);
  traceMetaData.latsta = srStation.latitude(ista);
  traceMetaData.lonevt = srEvent.longitude(ievt);
  traceMetaData.latevt = srEvent.latitude(ievt);
end
traceMetaData.nstack = ones(ntrace,1);
traceMetaData.sortIndex = (1:ntrace)';
traceMetaData.static = zeros(ntrace,1);
traceMetaData.mute = Inf(ntrace,1);


% Loop to read in traceDatas
fileIndexLast = -1;
for i = 1:ntrace
  if fileIndex(i)
    
    % Open new SEGY file

    if fileIndex(i)~=fileIndexLast
      j = fileIndex(i);
      fileIndexLast = j;
      if exist('fid')
        fclose(fid);
      end
      if segyCatalog.littleEndian(j)
          fid = fopen(cell2mat(segyCatalog.filename(j)),'r','ieee-le');
      else
          fid = fopen(cell2mat(segyCatalog.filename(j)),'r','ieee-be');
      end
      if segyCatalog.format(j)==1 || segyCatalog.format(j)==2 || segyCatalog.format(j)==5
        bytesPerSample = 4;
      elseif segyCatalog.format(j)==3  
        bytesPerSample = 2;
      elseif segyCatalog.format(j)==8  
        bytesPerSample = 1;
      else 
        error('get_segy.m - Reel format not valid/implemented') 
      end
    end
    
    % Get to start of traceMetaData
    if fseek(fid,byte0(i)+(i0(i)-1)*bytesPerSample,'bof')
      error('get_segy.m - Could not seek starting byte of trace')
    end
    
    % Read the traceMetaData
    if segyCatalog.format(j)==5
      traceData(1:nsamp(i),i) = fread(fid,nsamp(i),'float32');
    elseif segyCatalog.format(j)==2
      traceData(1:nsamp(i),i) = fread(fid,nsamp(i),'int32');
    elseif segyCatalog.format(j)==3
      traceData(1:nsamp(i),i) = fread(fid,nsamp(i),'int16');
    elseif segyCatalog.format(j)==8
      traceData(1:nsamp(i),i) = fread(fid,nsamp(i),'int8');
    elseif segyCatalog.format(j)==1
      c = fread(fid,nsamp(i)*4,'uchar');
      qf = ((c(4:4:nsamp(i)*4)/256+c(3:4:nsamp(i)*4))/256+c(2:4:nsamp(i)*4))/256;
      qc = c(1:4:nsamp(i)*4);
      qs = (qc<128)*2-1;
      qc = qc-128.*(qs<1)-64;
      traceData(1:nsamp(i),i) = qs.*16 .^(qc) .*qf;
    end
  end
end

% ETOMO Specific Start
% Create channel 5
vector = string2vector(paramSegy.channel);
if channel5
  scale = vector(find(vector==5)+1);
  scale = scale(end);
  i = true(size(vector));
  i(find(vector==5)+1) = false;
  vector = vector(i);
  index1 = find(traceMetaData.channel == 1);
  n = length(index1);
  traceMetaData.nsamp = [traceMetaData.nsamp; zeros(n,1)];
  traceMetaData.samprate = [traceMetaData.samprate; zeros(n,1)];
  traceMetaData.t0 = [traceMetaData.t0; zeros(n,1)];
  traceMetaData.eventid = [traceMetaData.eventid; zeros(n,1)];
  traceMetaData.eventIndex = [traceMetaData.eventIndex; zeros(n,1)];
  traceMetaData.stationIndex = [traceMetaData.stationIndex; zeros(n,1)];
  traceMetaData.channel = [traceMetaData.channel; zeros(n,1)];
  traceMetaData.range = [traceMetaData.range; zeros(n,1)];
  traceMetaData.azimuth = [traceMetaData.azimuth; zeros(n,1)];
  traceMetaData.xmid = [traceMetaData.xmid; zeros(n,1)];
  traceMetaData.ymid = [traceMetaData.ymid; zeros(n,1)];
  traceMetaData.xsta = [traceMetaData.xsta; zeros(n,1)];
  traceMetaData.ysta = [traceMetaData.ysta; zeros(n,1)];
  traceMetaData.xevt = [traceMetaData.xevt; zeros(n,1)];
  traceMetaData.yevt = [traceMetaData.yevt; zeros(n,1)];
  traceMetaData.lonmid = [traceMetaData.lonmid; zeros(n,1)];
  traceMetaData.latmid = [traceMetaData.latmid; zeros(n,1)];
  traceMetaData.lonsta = [traceMetaData.lonsta; zeros(n,1)];
  traceMetaData.latsta = [traceMetaData.latsta; zeros(n,1)];
  traceMetaData.lonevt = [traceMetaData.lonevt; zeros(n,1)];
  traceMetaData.latevt = [traceMetaData.latevt; zeros(n,1)];
  traceMetaData.nstack = [traceMetaData.nstack; zeros(n,1)];
  traceMetaData.sortIndex = [traceMetaData.sortIndex; zeros(n,1)];
  traceMetaData.static = [traceMetaData.static; zeros(n,1)];
  traceMetaData.mute = [traceMetaData.mute; zeros(n,1)];
  traceData = [traceData zeros(size(traceData,1),n)];
  k = traceMetaData.ntrace;
  for i = index1(:)'
    j = find(traceMetaData.channel(1:traceMetaData.ntrace) == 4 & ...
             traceMetaData.eventid(1:traceMetaData.ntrace) == traceMetaData.eventid(i) & ...
             strcmpi(traceMetaData.station(traceMetaData.ntrace),traceMetaData.station(i)));
    if isempty(j)
      error('Channel 1 & 4 not found for all station/eventid - this is required for creation of channel 5')
    end
    if ( traceMetaData.t0(i) ~= traceMetaData.t0(j) || ... 
         traceMetaData.samprate(i) ~= traceMetaData.samprate(j) || ...
         traceMetaData.nsamp(i) ~= traceMetaData.nsamp(j))
      error('Channel 1 & 4 sample or timing not matched - this is required for creation of channel 5')
    else
      k = k + 1;
      traceData(1:traceMetaData.nsamp(i),k) = (traceData(1:traceMetaData.nsamp(i),i) + ...
                                              scale * traceData(1:traceMetaData.nsamp(i),j)) / 2;
      traceMetaData.nsamp(k) = traceMetaData.nsamp(i);
      traceMetaData.samprate(k) = traceMetaData.samprate(i);
      traceMetaData.t0(k) = traceMetaData.t0(i);
      traceMetaData.eventid(k) = traceMetaData.eventid(i);
      traceMetaData.eventIndex(k) = traceMetaData.eventIndex(i);
      traceMetaData.station(k) = traceMetaData.station(i);
      traceMetaData.stationIndex(k) = traceMetaData.stationIndex(i);
      traceMetaData.channel(k) = 5;
      traceMetaData.range(k) = traceMetaData.range(i);
      traceMetaData.azimuth(k) = traceMetaData.azimuth(i);
      traceMetaData.xmid(k) = traceMetaData.xmid(i);
      traceMetaData.ymid(k) = traceMetaData.ymid(i);
      traceMetaData.xsta(k) = traceMetaData.xsta(i);
      traceMetaData.ysta(k) = traceMetaData.ysta(i);
      traceMetaData.xevt(k) = traceMetaData.xevt(i);
      traceMetaData.yevt(k) = traceMetaData.yevt(i);
      traceMetaData.lonmid(k) = traceMetaData.lonmid(i);
      traceMetaData.latmid(k) = traceMetaData.latmid(i);
      traceMetaData.lonsta(k) = traceMetaData.lonsta(i);
      traceMetaData.latsta(k) = traceMetaData.latsta(i);
      traceMetaData.lonevt(k) = traceMetaData.lonevt(i);
      traceMetaData.latevt(k) = traceMetaData.latevt(i);
      traceMetaData.nstack(k) = 2;      
      traceMetaData.sortIndex(k) = k;
      traceMetaData.static(k) = traceMetaData.static(i);
      traceMetaData.mute(k) = traceMetaData.mute(i);
    end
  end
  traceMetaData.ntrace = k;
end

% Create horizontal rotated
channel6or7 = false;
for ichan = [6 7 -6 -7];
  if any(vector==ichan)
    channel6or7 = true;
    index1 = find(traceMetaData.channel == 2);
    n = length(index1);
    traceMetaData.nsamp = [traceMetaData.nsamp; zeros(n,1)];
    traceMetaData.samprate = [traceMetaData.samprate; zeros(n,1)];
    traceMetaData.t0 = [traceMetaData.t0; zeros(n,1)];
    traceMetaData.eventid = [traceMetaData.eventid; zeros(n,1)];
    traceMetaData.eventIndex = [traceMetaData.eventIndex; zeros(n,1)];
    traceMetaData.stationIndex = [traceMetaData.stationIndex; zeros(n,1)];
    traceMetaData.channel = [traceMetaData.channel; zeros(n,1)];
    traceMetaData.range = [traceMetaData.range; zeros(n,1)];
    traceMetaData.azimuth = [traceMetaData.azimuth; zeros(n,1)];
    traceMetaData.xmid = [traceMetaData.xmid; zeros(n,1)];
    traceMetaData.ymid = [traceMetaData.ymid; zeros(n,1)];
    traceMetaData.xsta = [traceMetaData.xsta; zeros(n,1)];
    traceMetaData.ysta = [traceMetaData.ysta; zeros(n,1)];
    traceMetaData.xevt = [traceMetaData.xevt; zeros(n,1)];
    traceMetaData.yevt = [traceMetaData.yevt; zeros(n,1)];
    traceMetaData.lonmid = [traceMetaData.lonmid; zeros(n,1)];
    traceMetaData.latmid = [traceMetaData.latmid; zeros(n,1)];
    traceMetaData.lonsta = [traceMetaData.lonsta; zeros(n,1)];
    traceMetaData.latsta = [traceMetaData.latsta; zeros(n,1)];
    traceMetaData.lonevt = [traceMetaData.lonevt; zeros(n,1)];
    traceMetaData.latevt = [traceMetaData.latevt; zeros(n,1)];
    traceMetaData.nstack = [traceMetaData.nstack; zeros(n,1)];
    traceMetaData.sortIndex = [traceMetaData.sortIndex; zeros(n,1)];
    traceMetaData.static = [traceMetaData.static; zeros(n,1)];
    traceMetaData.mute = [traceMetaData.mute; zeros(n,1)];
    traceData = [traceData zeros(size(traceData,1),n)];
    k = traceMetaData.ntrace;
    for i = index1(:)'
      j = find(traceMetaData.channel(1:traceMetaData.ntrace) == 3 & ...
               traceMetaData.eventid(1:traceMetaData.ntrace) == traceMetaData.eventid(i) & ...
               strcmpi(traceMetaData.station(traceMetaData.ntrace),traceMetaData.station(i)));
      if isempty(j)
        error('Channel 2 & 3 not found for all station/eventid - this is required for creation of rotated horizontals')
      end
      if ( traceMetaData.t0(i) ~= traceMetaData.t0(j) || ... 
           traceMetaData.samprate(i) ~= traceMetaData.samprate(j) || ...
           traceMetaData.nsamp(i) ~= traceMetaData.nsamp(j))
        error('Channel 2 & 3 sample or timing not matched - this is required for creation of of rotated horizontals')
      end
      k = k + 1;
      l = find(strcmpi([srStation.name],traceMetaData.station(i)));
      if ichan>0
        rot2 = srStationRotation(l).rot;
        scale = srStationRotation(l).scale;
      else
        rot2 = srStationRotation(l).rotAlt;
        scale = srStationRotation(l).scaleAlt;
      end      
      if isempty(rot2) || isempty(scale)
        error('Horizontal Orientation data not available for all stations')
      end
      rot3 = rot2 - 90*sign(scale);
      direc2 = [sind(rot2); cosd(rot2)]; 
      direc3 = [sind(rot3); cosd(rot3)];
      if abs(ichan)==6;
        direcRad = [sind(traceMetaData.azimuth(i)); cosd(traceMetaData.azimuth(i))];
        scale2 = dot(direc2,direcRad) / sqrt(abs(scale));
        scale3 = dot(direc3,direcRad) * sqrt(abs(scale));
      elseif abs(ichan)==7;
        direcTran = [sind(traceMetaData.azimuth(i)-90); cosd(traceMetaData.azimuth(i)-90)];
        scale2 = dot(direc2,direcTran) / sqrt(abs(scale));
        scale3 = dot(direc3,direcTran) * sqrt(abs(scale)); 
      end
      traceData(1:traceMetaData.nsamp(i),k) = scale2 * traceData(1:traceMetaData.nsamp(i),i) + ...
                                              scale3 * traceData(1:traceMetaData.nsamp(i),j);
      traceMetaData.nsamp(k) = traceMetaData.nsamp(i);
      traceMetaData.samprate(k) = traceMetaData.samprate(i);
      traceMetaData.t0(k) = traceMetaData.t0(i);
      traceMetaData.eventid(k) = traceMetaData.eventid(i);
      traceMetaData.eventIndex(k) = traceMetaData.eventIndex(i);
      traceMetaData.station(k) = traceMetaData.station(i);
      traceMetaData.stationIndex(k) = traceMetaData.stationIndex(i);
      traceMetaData.channel(k) = ichan;
      traceMetaData.range(k) = traceMetaData.range(i);
      traceMetaData.azimuth(k) = traceMetaData.azimuth(i);
      traceMetaData.xmid(k) = traceMetaData.xmid(i);
      traceMetaData.ymid(k) = traceMetaData.ymid(i);
      traceMetaData.xsta(k) = traceMetaData.xsta(i);
      traceMetaData.ysta(k) = traceMetaData.ysta(i);
      traceMetaData.xevt(k) = traceMetaData.xevt(i);
      traceMetaData.yevt(k) = traceMetaData.yevt(i);
      traceMetaData.lonmid(k) = traceMetaData.lonmid(i);
      traceMetaData.latmid(k) = traceMetaData.latmid(i);
      traceMetaData.lonsta(k) = traceMetaData.lonsta(i);
      traceMetaData.latsta(k) = traceMetaData.latsta(i);
      traceMetaData.lonevt(k) = traceMetaData.lonevt(i);
      traceMetaData.latevt(k) = traceMetaData.latevt(i);
      traceMetaData.nstack(k) = traceMetaData.nstack(i);
      traceMetaData.sortIndex(k) = k;
      traceMetaData.static(k) = traceMetaData.static(i);
      traceMetaData.mute(k) = traceMetaData.mute(i);
    end
    traceMetaData.ntrace = k;
  end
end
  
% Discard channels 1-4 if only used for derived channels
if channel5 || channel6or7
  i = true(traceMetaData.ntrace,1);
  for ichan = 1:nchan
    if discardAtEnd(ichan)
      i(traceMetaData.channel==ichan) = false;
    end
  end
  traceData = traceData(:,i);
  traceMetaData.nsamp = traceMetaData.nsamp(i);
  traceMetaData.samprate = traceMetaData.samprate(i);
  traceMetaData.t0 = traceMetaData.t0(i);
  traceMetaData.eventid = traceMetaData.eventid(i);
  traceMetaData.eventIndex = traceMetaData.eventIndex(i);
  traceMetaData.station = traceMetaData.station(i);
  traceMetaData.stationIndex = traceMetaData.stationIndex(i);
  traceMetaData.channel = traceMetaData.channel(i);
  traceMetaData.range = traceMetaData.range(i);
  traceMetaData.azimuth = traceMetaData.azimuth(i);
  traceMetaData.xmid = traceMetaData.xmid(i);
  traceMetaData.ymid = traceMetaData.ymid(i);
  traceMetaData.xsta = traceMetaData.xsta(i);
  traceMetaData.ysta = traceMetaData.ysta(i);
  traceMetaData.xevt = traceMetaData.xevt(i);
  traceMetaData.yevt = traceMetaData.yevt(i);
  traceMetaData.lonmid = traceMetaData.lonmid(i);
  traceMetaData.latmid = traceMetaData.latmid(i);
  traceMetaData.lonsta = traceMetaData.lonsta(i);
  traceMetaData.latsta = traceMetaData.latsta(i);
  traceMetaData.lonevt = traceMetaData.lonevt(i);
  traceMetaData.latevt = traceMetaData.latevt(i);
  traceMetaData.ntrace = sum(i);
  traceMetaData.sortIndex = (1:sum(i))';
  traceMetaData.static = traceMetaData.static(i);
  traceMetaData.mute = traceMetaData.mute(i);
end

% Create Stacked Undershoot Lines
vector = string2vector(paramSegy.eventid);
if any(vector == -52) || any(vector == -55) || any(vector == -53)
  if reductionSlowness
    error('Events -52, -53 and -55: Cannot load stacked ETOMO SEGY data with a reduction velocity');
  end
  keep = true(1,traceMetaData.ntrace);
  if any(vector == -52)  
    for i = 18:259
      for channel = 1:5
        for station = traceMetaData.stationList
          j = find((traceMetaData.eventid==7000+i | traceMetaData.eventid==8000+i | traceMetaData.eventid==39000+i) & ...
            traceMetaData.channel==channel & strcmp(traceMetaData.station,station));
          if length(j) > 1
            if any(diff(traceMetaData.nsamp(j))) || any(diff(traceMetaData.samprate(j)))
              error('Trace mismatch for stacked line -52')
            else
              if any(diff(traceMetaData.t0(j)))
                disp(['Trace start time offset of ' num2str(diff(traceMetaData.t0(j))') ' for shots ' int2str(traceMetaData.eventid(j)')])
              end
              traceData(:,j(1)) = sum(traceData(:,j)')' / length(j);
              %traceMetaData.eventid(j(1)) = 52000 + i;
              traceMetaData.lonevt(j(1)) = mean(traceMetaData.lonevt(j));
              traceMetaData.latevt(j(1)) = mean(traceMetaData.latevt(j));
              traceMetaData.nstack(j(1)) = traceMetaData.nstack(j(1))*length(j);
              keep(j(2:end)) = false;
              disp(['Stacked traces ' int2str(traceMetaData.eventid(j)')])
            end
          end
        end
      end
    end
  end
  if any(vector == -55)  
% Fix offset in shot numbers (i.e., line 1 shot 1021 = line 42 shot 42020) 
%     for i = 13:280
    for i = 14:280
      for channel = 1:5
        for station = traceMetaData.stationList
%           j = find((traceMetaData.eventid==1000+i | traceMetaData.eventid==42000+i) & ...
%             traceMetaData.channel==channel & strcmp(traceMetaData.station,station));
          j = find((traceMetaData.eventid==1000+i | traceMetaData.eventid==42000+i-1) & ...
            traceMetaData.channel==channel & strcmp(traceMetaData.station,station));
          if length(j) > 1
            if any(diff(traceMetaData.nsamp(j))) || any(diff(traceMetaData.samprate(j)))
              error('Trace mismatch for stacked line -55')
            else
              if any(diff(traceMetaData.t0(j)))
                disp(['Trace start time offset of ' num2str(diff(traceMetaData.t0(j))') ' for shots ' int2str(traceMetaData.eventid(j)')])
              end
              traceData(:,j(1)) = sum(traceData(:,j)')' / length(j);
              %traceMetaData.eventid(j(1)) = 55000 + i;
              traceMetaData.lonevt(j(1)) = mean(traceMetaData.lonevt(j));
              traceMetaData.latevt(j(1)) = mean(traceMetaData.latevt(j));
              traceMetaData.nstack(j(1)) = traceMetaData.nstack(j(1))*length(j);
              keep(j(2:end)) = false;
              disp(['Stacked traces ' int2str(traceMetaData.eventid(j)')])
            end
          end
        end
      end
    end
  end
  if any(vector == -53)  
    for i = 21:264
      for channel = 1:5
        for station = traceMetaData.stationList
          j = find((traceMetaData.eventid==4000+i | traceMetaData.eventid==6000+i) & ...
            traceMetaData.channel==channel & strcmp(traceMetaData.station,station));
          if length(j) > 1
            if any(diff(traceMetaData.nsamp(j))) || any(diff(traceMetaData.samprate(j)))
              error('Trace mismatch for stacked line -53')
            else
              if any(diff(traceMetaData.t0(j)))
                disp(['Trace start time offset of ' num2str(diff(traceMetaData.t0(j))') ' for shots ' int2str(traceMetaData.eventid(j)')])
              end
              traceData(:,j(1)) = sum(traceData(:,j)')' / length(j);
              %traceMetaData.eventid(j(1)) = 53000 + i;
              traceMetaData.lonevt(j(1)) = mean(traceMetaData.lonevt(j));
              traceMetaData.latevt(j(1)) = mean(traceMetaData.latevt(j));
              traceMetaData.nstack(j(1)) = traceMetaData.nstack(j(1))*length(j);
              keep(j(2:end)) = false;
              disp(['Stacked traces ' int2str(traceMetaData.eventid(j)')])
            end
          end
        end
      end
    end
  end
  traceData = traceData(:,keep);
  traceMetaData.nsamp = traceMetaData.nsamp(keep);
  traceMetaData.samprate = traceMetaData.samprate(keep);
  traceMetaData.t0 = traceMetaData.t0(keep);
  traceMetaData.eventid = traceMetaData.eventid(keep);
  traceMetaData.eventIndex = traceMetaData.eventIndex(keep);
  traceMetaData.station = traceMetaData.station(keep);
  traceMetaData.stationIndex = traceMetaData.stationIndex(keep);
  traceMetaData.channel = traceMetaData.channel(keep);
  traceMetaData.range = traceMetaData.range(keep);
  traceMetaData.azimuth = traceMetaData.azimuth(keep);
  traceMetaData.xmid = traceMetaData.xmid(keep);
  traceMetaData.ymid = traceMetaData.ymid(keep);
  traceMetaData.xsta = traceMetaData.xsta(keep);
  traceMetaData.ysta = traceMetaData.ysta(keep);
  traceMetaData.xevt = traceMetaData.xevt(keep);
  traceMetaData.yevt = traceMetaData.yevt(keep);
  traceMetaData.lonmid = traceMetaData.lonmid(keep);
  traceMetaData.latmid = traceMetaData.latmid(keep);
  traceMetaData.lonsta = traceMetaData.lonsta(keep);
  traceMetaData.latsta = traceMetaData.latsta(keep);
  traceMetaData.lonevt = traceMetaData.lonevt(keep);
  traceMetaData.latevt = traceMetaData.latevt(keep);
  traceMetaData.nstack = traceMetaData.nstack(keep);
  traceMetaData.sortIndex = (1:sum(keep))';
  traceMetaData.static = traceMetaData.static(keep);
  traceMetaData.mute = traceMetaData.mute(keep);
  traceMetaData.ntrace = sum(keep);
end
% ETOMO Specific End

if exist('fid')
  fclose(fid);
end  
