% Script to orient horizontal geophones 
%   This script is included with tlPicker but not accessible from the GUI menu
%
%   The script orients the horizontals using a brute force grid search and 
%     searches for an orientation that minimizes the mean misfit of
%     particle motions assuming both uniform scaling on the two horizontal
%     channels and differential scaling.  If the vertical channel is
%     available it attempts to use motions for shots at short ranges to
%     eliminate a 180° ambiguity in the direction of positive motion.
%  
%  Parameters for the orientation are set within a field orientHoriz of
%    currentMenu (or defaultMenu).  It has the following fields
%    srStationRotation - File in which to store results and from which
%                        they are accessed to create radial and tangential
%                        horizontal channels 
%    tlim0    -  Start of window relative to water pick for determinining 
%                particle motions using PCA (normally 0 s)
%    tlim1    -  End of window relative to water pick for determinining 
%                particle motions using PCA (0.05 s used for ETOMO)
%    minEvt   -  Minimum number of events required to orient horizontals
%    minScale  
%    maxScale
%    nScale   -  When determining the particle motion with differential
%                scaling of the horizontals the algorithm will consider
%                nScale different scalings of channel 3 logarthimically
%                spaced between minScale and maxScale.
%    rangeCrit - Range beneath which three component motions can reliably
%                give direction of arrival (it is basically the critical
%                range) (km)
%    plotEachEvt - Logical to determine if plots are made of each particle
%                  motion analysis
%    noVertical  - logical indicates that no vertial data is available
%                  If this is the case there is a 180° ambiguity in the
%                  orientation
%
%   Before running orientHoriz the data should loaded and plotted with the
%   following parameter settings
%      One station
%      All shots
%      Channels [1:3] (or [1:2] if noVertical == true)
%      Time limits [0 5]
%      Reduction velocity of 1.5
%      Range selection [0 8]
%      Selection DFS/Phase for water wave
%      Filter limits [40 80] works well
%      Static Pick DFS/Phase for water wave
%      Active Pick DFS/Phase for water wave
%
%  File srStationRotation holds one structure srStationRotation with one
%  element per station ordered as in srStation and the following fields
%    name      - Station name
%    use       - Logical to indicate if rotation is reliable
%    rot       - Azimuth of channel 2 relative to north assuming both
%                 horizontal channels have the same scaling
%    scale     - Scaling of channel 3.  If it is 1.0 then channel 3 is 90
%                 degrees anticlockwise from channel 2 looking down on
%                 instrument (i.e., a right hand coordinate system).  If it
%                 is -1.0 then channel 3 is clockwise of channel 2
%  The following fields are optional for tlPicker in that they are not
%  needed to rotate horizontals into radial and tangential channels
%    processed - Logical to indicate that this station has been processed
%                 for horizontal channel orientations - Useful to indicate
%                 an unsuccessful attempt to orient channes
%    meanMisfit - Mean absolute misfit of rotated particle motions (degree)
%    medianMisfit - Median absolute misfit of rotated particle motions (degree)
%    rotAlt     - Alternative azimuth of channel 2 based on allowing
%                 differential scaling
%    scaleAlt   - Alternative scaling for channel3 with sign convention
%                 as for scale 
%    meanMisfitAlt - Mean absolute misfit
%    medianMisfitAlt - mean absolute misfit;
%    segy        - currentMenu.segy;
%    filter      - currentMenu.filter;
%    static      - currentMenu.static;
%    active      - currentMenu.active;
%    orientHoriz - currentMenu.orientHoriz;

% Check appropriate data is loaded 
status = 0;
if length(traceMetaData.stationList)>1
  status=1;
  warning('orientHoriz aborting because data is loaded for more than one station')
  return
end
if isempty(traceMetaData.stationList)
  status=2;
  warning('orientHoriz aborting because no data is loaded')
  return
end
if active.npick~=traceMetaData.ntrace
  status=3;
  warning('orientHoriz aborting because pick file does not match trace metadata')
  return
end
if ~strcmpi(deblank_fb(currentMenu.static.DFS),deblank_fb(currentMenu.active.directory)) || ...
   ~strcmpi(deblank_fb(currentMenu.static.phase),deblank_fb(currentMenu.active.phase))
  status = 4;
  warning('orientHoriz aborting because active picks are not used for statics')
  return
end

% Event ids
eventids = sort(unique(traceMetaData.eventid));

j = 0;
clear rangePM baz c2 c3 xPM yPM

% Cycle through eventids
for eventid = eventids(:)';
  
  % Get data for 3 channels
  if ~currentMenu.orientHoriz.noVertical
    i1 = find(traceMetaData.eventid == eventid & traceMetaData.channel==1);
  end
  i2 = find(traceMetaData.eventid == eventid & traceMetaData.channel==2);
  i3 = find(traceMetaData.eventid == eventid & traceMetaData.channel==3);
  
  % Check data is present for each channel
  if ~currentMenu.orientHoriz.noVertical
    if length(i1)~=1 
      warning(['orientHoriz for eventid = ' int2str(eventid) ' has no data for channel 1 - Skipping'])
      continue
    end
  end
  if length(i2)~=1 
    warning(['orientHoriz for eventid = ' int2str(eventid) ' has no data for channel 2 - Skipping'])
    continue
  elseif length(i3)~=1 
    warning(['orientHoriz for eventid = ' int2str(eventid) ' has no data for channel 3 - Skipping'])
    continue
  end
  
  % Verify data is consistent
  if currentMenu.orientHoriz.noVertical
    i = [i2 i3];
  else
    i = [i1 i2 i3];
  end
  if any(diff(traceMetaData.t0(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent start time - Skipping'])
    continue
  end
  if any(diff(traceMetaData.samprate(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent sample rate - Skipping'])
    continue
  end
  if any(diff(traceMetaData.nsamp(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent number of samples - Skipping'])
    continue
  end
  if any(diff(traceMetaData.range(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent range - Skipping'])
    continue
  end
  if any(diff(traceMetaData.azimuth(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent azimuth - Skipping'])
    continue
  end  
  if any(isnan(active.time(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - No Picks - Skipping'])
    continue
  end
  if any(diff(active.time(i)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Channels have inconsistent picks - Skipping'])
    continue
  end
     
  index0 = round((active.time(i2)-traceMetaData.t0(i2)+currentMenu.orientHoriz.tlim0)*traceMetaData.samprate(i2) + 1);
  if index0<1 || index0>traceMetaData.nsamp(i2);
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Start of particle time window not within loaded trace - Skipping'])
    continue
  end
  index1 = round((active.time(i2)-traceMetaData.t0(i2)+currentMenu.orientHoriz.tlim1)*traceMetaData.samprate(i2) + 1);
  if index1<1 || index1>traceMetaData.nsamp(i2);
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - End of particle time window not within loaded trace - Skipping'])
    continue
  end
  data = traceData(index0:index1,i);
  if currentMenu.orientHoriz.noVertical
    data = [zeros(size(data,1),1) data];
  end  
  if any(any(isnan(data)))
    warning(['orientHoriz for eventid = ' int2str(eventid) ' - Not all requested time window is loaded - Skipping'])
    continue
  end
  
  for k=1:3
    data(:,k) = data(:,k) - mean(data(:,k));
  end
  
  % PCA of horizontal motions 
  % 2D Horizontal
  [V2,D2] = eig(data(:,2:3)'*data(:,2:3));
  D2 = diag(D2);
  [~,k] = max(D2);
  frac2 = D2(k)/sum(D2);
  V2 = V2(:,k);
  scale = max(max(abs(data(:,2:3))));
  
  % 3D Horizontal
  if ~currentMenu.orientHoriz.noVertical
    [V3,D3] = eig(data'*data);
    D3 = diag(D3);
    [~,k] = max(D3);
    frac3 = D3(k)/sum(D3);
    V3 = V3(:,k);
  
    % Adjust V2 so it is the back azimuth for a downward motion -
    % Will only work reliably for incidence angles below critical angle
    if V3(1)>0
      if dot(V3(2:3),V2(1:2))<0
        V2 = - V2;
      end
    elseif V3(1)<0
      if dot(V3(2:3),V2(1:2))>0
        V2 = - V2;
      end
    end  
  end
  
  % Get back azimuth for this event
  j = j+1;
  baz(j) = traceMetaData.azimuth(i2)-180;
  if baz(j)<0; baz(j) = baz(j) + 360; end
  
  % Components of water wave particle motion in directions of two
  % horizontal channels and x/y/range values of event
  c2(j) = V2(1);
  c3(j) = V2(2);
  xPM(j) = traceMetaData.xevt(i2);
  yPM(j) = traceMetaData.yevt(i2);
  rangePM(j) = traceMetaData.range(i2);

%% Plotting each determination
  if currentMenu.orientHoriz.plotEachEvt
    %Horizontal PM
    figure(103)
    clf
    plot(data(:,2),data(:,3),data(1,2),data(1,3),'bo','markersize',10,'markerfacecolor','b');
    axis equal
    hold on
    plot([V2(1) -V2(1)]*scale,[V2(2) -V2(2)]*scale,'r')
    if ~currentMenu.orientHoriz.noVertical
      plot([V3(2) -V3(2)]*scale,[V3(3) -V3(3)]*scale,'g')
    end
    xlabel('X, km')
    ylabel('Y, km')
    title(['Waveforms - eventid=' int2str(traceMetaData.eventid(i2)) '  Range=' num2str(traceMetaData.range(i2)) ...
            '  Azimuth=' num2str(traceMetaData.azimuth(i2)) ]);

    % Waveforms
    figure(104)
    clf
    t = currentMenu.orientHoriz.tlim0+(0:size(data,1)-1)/traceMetaData.samprate(i1);
    plot(t,data)
    xlabel('Time, s')
    ylabel('Counts')
    title(['Waveforms - eventid=' int2str(traceMetaData.eventid(i2)) '  Range=' num2str(traceMetaData.range(i2)) ...
            '  Azimuth=' num2str(traceMetaData.azimuth(i2)) ]);
    pause 
  end

end

%% Plot Apparent back azimuth map
figure(105); clf;
plot(traceMetaData.xsta(1),traceMetaData.ysta(1),'ko','markersize',16,'markerfacecolor','k')
xlabel('X, km')
ylabel('Y, km')
title(['Observations Station=' cell2mat(traceMetaData.stationList)]);
axis equal
hold on
for j = 1:length(xPM)
  plot(xPM(j)+[0 -c3(j)],yPM(j)+[0 c2(j)],'-r',xPM(j),yPM(j),'or');
%   text(xPM(j),yPM(j),[int2str(round(baz(j)))]);
end

if length(xPM)<currentMenu.orientHoriz.minEvt
  warning('orientHoriz aborting because not enough particle motion data')
  return 
end

%% Find Rotation assuming scale is uniform on horizontal channels
% Consider channel 3 on either side of channel 2.
% For a positive scaling of channel 2 the assumption is that channel 3 is
% 90 degrees anticlockwise from channel 2 looking down on OBS
meanMisfitFixedScale = inf;
for rotTry = 0:1:179;
  for scaleTry = [-1:2:1]
    bazPred = rotTry + atan2(-c3*scaleTry,c2)*180/pi;
    e = baz - bazPred;
    e(e>90) = e(e>90)-180;
    e(e>90) = e(e>90)-180;
    e(e<-90) = e(e<-90)+180;
    e(e<-90) = e(e<-90)+180;
    meanMisfit = mean(abs(e));
    if meanMisfit < meanMisfitFixedScale
      rotFixedScale = rotTry;
      scaleFixedScale = scaleTry;
      medianMisfitFixedScale = median(abs(e));
      meanMisfitFixedScale = meanMisfit;
    end
  end
end

% Adjust by 180 degrees to give direction of positive motion
if ~currentMenu.orientHoriz.noVertical
  bazPred = rotFixedScale + atan2(-c3*scaleFixedScale,c2)*180/pi;
  misfit = [bazPred(rangePM<currentMenu.orientHoriz.rangeCrit)' - baz(rangePM<currentMenu.orientHoriz.rangeCrit)'];
  misfit(misfit>180) = misfit(misfit>180) - 360;
  misfit(misfit<-180) = misfit(misfit<-180) + 360;
  if median(abs(misfit))>90
    rotFixedScale = rotFixedScale+180;
    misfit = misfit + 180;
    misfit(misfit>270) = misfit(misfit>270)-360;
  end
  disp(['Fixed scaling - 180 degree ambiguity eliminated with ' int2str(sum(abs(misfit)<90)) ' of ' int2str(length(misfit)) ' short range observations'])
  ambiguityFixed = [sum(abs(misfit)<90) length(misfit)];
else
  disp('Fixed scaling - 180 degree ambiguity not eliminated')
  ambiguityFixed = [0 0];
end

% Plot uniform scaling results predictions
figure(106); clf;
plot(traceMetaData.xsta(1),traceMetaData.ysta(1),'ko','markersize',16,'markerfacecolor','k')
xlabel('X, km')
ylabel('Y, km')
title(['Fixed Scaling Model Station=' cell2mat(traceMetaData.stationList) '  Rotation=' int2str(rotFixedScale) '  Scale=' int2str(scaleFixedScale) '  Mean Misfit=' num2str(meanMisfitFixedScale)]);
axis equal
hold on
for deg=0:15:360
  plot(traceMetaData.xsta(1)+[0 sind(deg)]*10,traceMetaData.ysta(1)+[0 cosd(deg)]*10,'k-');
end
bazPred = rotFixedScale + atan2(-c3*scaleFixedScale,c2)*180/pi;
bazPred(bazPred-baz>90) = bazPred(bazPred-baz>90)-180;
bazPred(bazPred-baz>90) = bazPred(bazPred-baz>90)-180;
bazPred(bazPred-baz<-90) = bazPred(bazPred-baz<-90)+180;
bazPred(bazPred-baz<-90) = bazPred(bazPred-baz<-90)+180;
bazPlot = bazPred + srGeometry.rotation;
for j = 1:length(xPM);
  plot(xPM(j)+[0 sind(bazPlot(j))],yPM(j)+[0 cosd(bazPlot(j))],'b-')
end

%% Find Rotation assuming scale is not uniform on horizontal channels
meanMisfitFreeScale = inf;
for rotTry = 0:1:179;
  for scaleTry = logspace(log10(currentMenu.orientHoriz.minScale),log10(currentMenu.orientHoriz.maxScale),currentMenu.orientHoriz.nScale)*sign(scaleFixedScale)
    bazPred = rotTry + atan2(-c3*scaleTry,c2)*180/pi;
    e = baz - bazPred;
    e(e>90) = e(e>90)-180;
    e(e>90) = e(e>90)-180;
    e(e<-90) = e(e<-90)+180;
    e(e<-90) = e(e<-90)+180;
    meanMisfit = mean(abs(e));
    if meanMisfit < meanMisfitFreeScale
      rotFreeScale = rotTry;
      scaleFreeScale = scaleTry;
      medianMisfitFreeScale = median(abs(e));
      meanMisfitFreeScale = meanMisfit;
    end
  end
end

% Adjust by 180 degrees to give direction of positive motion
if ~currentMenu.orientHoriz.noVertical
  bazPred = rotFreeScale + atan2(-c3*scaleFreeScale,c2)*180/pi;
  misfit = [bazPred(rangePM<currentMenu.orientHoriz.rangeCrit)' - baz(rangePM<currentMenu.orientHoriz.rangeCrit)'];
  misfit(misfit>180) = misfit(misfit>180) - 360;
  misfit(misfit<-180) = misfit(misfit<-180) + 360;
  if median(abs(misfit))>90
    rotFreeScale = rotFreeScale+180;
    misfit = misfit + 180;
    misfit(misfit>180) = misfit(misfit>180)-360;
  end
  disp(['Free scaling - 180 degree ambiguity eliminated with ' int2str(sum(abs(misfit)<90)) ' of ' int2str(length(misfit)) ' short range observations'])
  ambiguityFree = [sum(abs(misfit)<90) length(misfit)];
else
  disp('Free scaling - 180 degree ambiguity not eliminated')
  ambiguityFree = [0 0];
end


% Plot non uniform scaling results predictions
figure(107); clf;
plot(traceMetaData.xsta(1),traceMetaData.ysta(1),'ko','markersize',16,'markerfacecolor','k')
xlabel('X, km')
ylabel('Y, km')
title(['Free Scaling Model Station=' cell2mat(traceMetaData.stationList) '  Rotation=' int2str(rotFreeScale) '  Scale=' num2str(scaleFreeScale) '  Mean Misfit=' num2str(meanMisfitFreeScale)]);
axis equal
hold on
for deg=0:15:360
  plot(traceMetaData.xsta(1)+[0 sind(deg)]*10,traceMetaData.ysta(1)+[0 cosd(deg)]*10,'k-');
end
bazPred = rotFreeScale + atan2(-c3*scaleFreeScale,c2)*180/pi;
bazPred(bazPred-baz>90) = bazPred(bazPred-baz>90)-180;
bazPred(bazPred-baz>90) = bazPred(bazPred-baz>90)-180;
bazPred(bazPred-baz<-90) = bazPred(bazPred-baz<-90)+180;
bazPred(bazPred-baz<-90) = bazPred(bazPred-baz<-90)+180;
bazPlot = bazPred + srGeometry.rotation;
for j = 1:length(xPM);
  plot(xPM(j)+[0 sind(bazPlot(j))],yPM(j)+[0 cosd(bazPlot(j))],'b-')
end 

%% Save rotation
i = find(strcmpi([srStationRotation.name], traceMetaData.stationList));
if ~srStationRotation(i).processed
  disp('No rotation information is currently saved for this station')
  proceed = input('Do you want to save this new rotation (y/N) : ','s'); 
else
  if ~srStationRotation(i).use
    disp('Rotation information has previously been saved for this station but is flagged as unusable');
  else
    disp('Rotation information has previously been saved for this station');
  end
  proceed = input('Do you want to overwrite with this new rotation (y/N) : ','s');
end  
proceed = strcmpi(proceed,'y');
if proceed
  use = input('Do you want to flag this new rotation as useable (y/N) : ','s');
  use = strcmpi(use,'y');
  i = find(strcmpi([srStationRotation.name], traceMetaData.stationList));
  srStationRotation(i).name = srStation.name(i);
  srStationRotation(i).use = use;
  srStationRotation(i).rot = rotFixedScale;
  srStationRotation(i).scale = scaleFixedScale;
  % Optional
  srStationRotation(i).ambiguity = ambiguityFixed;
  srStationRotation(i).processed = true;
  srStationRotation(i).meanMisfit = meanMisfitFixedScale;
  srStationRotation(i).medianMisfit = medianMisfitFixedScale;
  srStationRotation(i).rotAlt = rotFreeScale ;
  srStationRotation(i).scaleAlt = scaleFreeScale;
  srStationRotation(i).ambiguityAlt = ambiguityFree;
  srStationRotation(i).meanMisfitAlt = meanMisfitFreeScale;
  srStationRotation(i).medianMisfitAlt = medianMisfitFreeScale;
  srStationRotation(i).segy = currentMenu.segy;
  srStationRotation(i).filter = currentMenu.filter;
  srStationRotation(i).static = currentMenu.static;
  srStationRotation(i).active = currentMenu.active;
  srStationRotation(i).orientHoriz = currentMenu.orientHoriz;
  [~,srStationRotation(i).user] = unix('whoami');
  eval(['save ' currentMenu.orientHoriz.srStationRotation ' srStationRotation'])
end
