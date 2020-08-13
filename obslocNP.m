function s = obslocNP(p)
% Function to relocate OBSs (and shots)
%prevents all printing of figures when doing the inversion 
%created by Maleen Kidiwela 05/05/2020

% Usage
%   s = obsloc(p)
%
% Inputs
%   p is a structure with the following fields
%    srControl   - Stingray control structure (only field tf_latlon is relevant)
%    srGeometry  - Stingray geometry structure
%    srStation   - Stringray station structure with fields x,y,z populated
%    srEvent     - Stingray event structure with fields x,y populated
%    tlArrival   - Stingray/TomoLab arrival structure with Pw phases;
%    srElevation - Stingray elevation structure for bathymetry
%    solve       - Structure controlling what is solved for with the following fields
%      xyStation       - solve for station x & y
%      zStation        - solve for station z
%      xyEvent         - solve for shot x & y
%      tEvent          - solve for shot instant adjustment
%      useBathymetry   - adjust apriori Station depths to bathymetry at
%                          current locations
%      adjustDataError - Assume travel time errors are relative and adjust absolute values to 
%                          match observed misfits when calculating errors
%    Station     - Structure containing the following fields
%      xErrorApr     - Scaler apriori station x error (km)
%      yErrorApr     - Scaler apriori station y error (km)
%      zErrorApr     - Scaler apriori station z error (km).  
%        Uses initial locations as apriori values and employs a jumping strategy
%    Event       - Structure containing the following fields
%      xErrorApr     - Scaler apriori shot x error (km)
%      yErrorApr     - Scaler apriori shot y error (km)
%      tErrorApr     - Scaler apriori shot instant error (km)
%        Since shot locations/instants can be badly constrained when OBS spacing is large 
%        these provide important damping when shots are included in inversion 
%    nlookup     - Number of files with travel time and derivative lookup tables
%    lookup           - Structure with nlookup elements and fields
%      file  - File name for a travel time / derivative lookup table file with the 
%              following variables (TTtables_simple.m creates file for constant velocity)
%        zTable    - Vector of Z (depth) values (km)
%        rTable    - Vector of horizontal range values (km)
%        ttTable   - Matrix of travel times (first index is range) (s)
%        dtdrTable - Matrix of horizontal travel time derivatives (s/km)
%        dtdzTable - Matrix of vertical travel time derivatives (s/km)
%        ztable    - Vector of Z (depth) values (km)
%      eventElevation - Shot elevation this lookup table is constructed for 
%    nlsqr       - Use LSQR when there are this many model parameters
%    iterMax    - Maximum number of inversion terations      
%    dtMax       - Convergence requires all shot instant changes are less than this     
%    dxyzMax     - Convergence requires all distance changes are less than this 
%    drmsNax     - Convergence requires change in the RMS to be less than
%                    this for two iterations.
%    drmsNormMax - Convergence requires change in the normalized RMS to be 
%                    less than this for two iterations
%    rangeMax    - Maximum range to use in the inversion
%    nStaPlot    - Number of stations on residual versus event index plots
%                  (0 for no plots)
%    plotLabel   - Logical to control addition of shot/station labels to plots
%    yLabelRMms  - Vertical position of evant/station labels (OPTIONAL)
%    yLabelRmsNorm  - Vertical position of event/station labels (OPTIONAL)
%
% Outputs
%   s is a structure with the following fields
%    xStation       - Relocated station x coordinates (km)
%    yStation       - Relocated station y coordinates (km)
%    zStation       - Relocated station z coordinates (km)
%    xStation_int      - Initial station x coordinates (km)
%    yStation_int      - Initial station y coordinates (km)
%    zStation_int      - Initial station z coordinates (km)
%    xEvent         - Relocated shot x coordinates (km)
%    yEvent         - Relocated shot y coordinates (km)
%    tEvent         - Relocated shot instant correction (s)
%    xStationApr    - Apriori station x coordinates (km)
%    yStationApr    - Apriori station y coordinates (km)
%    zStationApr    - Apriori station z coordinates (km)
%                     if p.useBathymetry==1 this will be the bathymetric
%                     depth used for the final iteration
%    xEventApr      - Apriori shot x coordinates (km)
%    yEventApr      - Apriori shot y coordinates (km)
%    tEventApr      - Apriori shot instant correction (s)
%    xStationError  - Uncertainty in relocated station x coordinates (km)
%    yStationError  - Uncertainty in relocated station y coordinates (km)
%    zStationError  - Uncertainty in relocated station y coordinates (km)
%    xEventError    - Uncertainty in relocated shot x coordinates (km)
%    yEventError    - Uncertainty in relocated shot y coordinates (km)
%    tEventError    - Uncertainty in relocated shot instant corrections (s)
%      All errors are set to zero if not included in inversion or using LSQR
%    rms            - RMS travel time misfit
%    rmsNorm        - RMS of travel time misfits normalized to travel time errors
%    rmsEvent       - Vector of RMS travel time misfits for each shot
%    rmsEventNorm   - Vector of RMS normalized travel time misfits for each shot
%    rmsStation     - Vector of RMS travel time misfits for each station
%    rmsStationNorm - Vector of RMS normalized travel time misfits for each station
%    ista           - vector of indicies in input srStation structure of stations included in inversion 
%    ievt           - vector of indicies in input srEvent structure of shots included in inversion
%    time           - Matrix of observed travel times with the first index representing
%                     the shot (ievt) and the second index representing the Station
%                    (ista).  NaNs indicate no data.  Data obtained from tlArrival structure
%    timePred       - Matrix of predicted travel times (s).  NaNs indicate no prediction
%    timeError      - Matrix of travel time errors (s).  Data obtained from tlArrival structure
%    srStation      - Stringray station structure with relocated stations
%    srEvent        - Stringray station structure with relocated shots
%
% The inversion generates a series of plots with figure numbers as follows
%    11       - X-Y plan view of stations and events
%    12,13    - X-Z and Y-Z cross sections (not very useful)
%    14       - Z versus station index
%    21       - Inversion RMS residual progression
%    22       - RMS residuals versus event index
%    23       - RMS residual versus station index
%    31-##    - Residuals versus event index (6/12 stations per plot)
%    51-##    - Normalized residuals versus event index (6/12 stations per plot)
%    71-73    - Change in station parameters for final interation
%    74-76    - Change in event parameters for final interation
%    77-79    - Change in station parameters for final interation from original
%
% William Wilcock & Dax Soule, February 2010 


%% Travel time and travel time uncertainty matrices
time = NaN(p.srEvent.nevt,p.srStation.nsta);
timeError = NaN(p.srEvent.nevt,p.srStation.nsta);
for i = 1:p.srStation.nsta
  sA = subset_tlArrival(p.tlArrival, 'Pw', char(p.srStation.name(i)), [], []);
  if ~isempty(sA.time);
    [j,status] = index_eventid(p.srEvent.id,sA.eventid);
    if status
      error(['index_eventid; status = ' int2str(status)])
    end
    time(j,i) = sA.time;
    timeError(j,i) = sA.error;
  end
end

% Select stations and events with data
inosta = find(~sum(~isnan(time)));
inoevt = find(~sum(~isnan(time')));
ista = find(sum(~isnan(time)))';
ievt = find(sum(~isnan(time')))';
%ista = ista([1:41,43:50,52:length(ista)],:); % Added since there are issues with OBSs 45 & 54
time = time(ievt,ista);
timeError = timeError(ievt,ista);

% Print warnings of missing stations and events
if ~isempty(inosta)
  fprintf('No data for Stations\n')
  for i=1:length(inosta)
    fprintf(['  ',cell2mat(p.srStation.name(inosta(i))),'\n'])
  end
  fprintf('\n')
end
if ~isempty(inosta)
  fprintf('No data for Events\n')
  for i=1:length(inoevt)
    fprintf('  %i\n',p.srEvent.id(inoevt(i)))
  end
  fprintf('\n')
end 

% Number of Stations and shots and indexing vectors
m = length(ievt);
n = length(ista);
[shotIndex,StationIndex] = find(ones(m,n));

%% Starting values for model parameters and apriori uncertainties
xStation = p.srStation.x(ista);
yStation = p.srStation.y(ista);
zStation = -p.srStation.elevation(ista);
zStation0 = -p.srStation.elevation(ista);
xEvent = p.srEvent.x(ievt);
yEvent = p.srEvent.y(ievt);
tEvent = zeros(length(ievt),1);
xStation_int = xStation;
yStation_int = yStation;
zStation_int = zStation;
xEvent_int = xEvent;
yEvent_int = yEvent;
tEvent_int = tEvent;

xStationErrorApr = zeros(n,1) + p.Station.xErrorApr;
yStationErrorApr = zeros(n,1) + p.Station.yErrorApr;
zStationErrorApr = zeros(n,1) + p.Station.zErrorApr;
xEventErrorApr   = zeros(m,1) + p.Event.xErrorApr;
yEventErrorApr   = zeros(m,1) + p.Event.yErrorApr;
tEventErrorApr   = zeros(m,1) + p.Event.tErrorApr;

%% Load travel time & spatial derivative lookup tables
lookup = ones(length(ievt),1);
for i = 1:p.nlookup
  eval(['load ' p.lookup(i).file])
  p.lookup(i).ttTable   = ttTable;
  p.lookup(i).dtdrTable = dtdrTable;
  p.lookup(i).dtdzTable = dtdzTable;
  p.lookup(i).rTable    = rTable;
  p.lookup(i).zTable    = zTable;
  if i>1
    lookup(p.srEvent.elevation(ievt)==p.lookup(i).eventElevation) = i;
  end
end

%% Iteration loop
converged = 0;
for iter = 1:p.iterMax+1;  
    
  % Forward problem with respect to new model parameters 
  DX = xEvent(shotIndex,:) - xStation(StationIndex,:);
  DY = yEvent(shotIndex,:) - yStation(StationIndex,:);
  range = sqrt( DX.^2 + DY.^2 );
  depth = zStation(StationIndex,:);
  timePred = NaN(m,n);
  dtdr = NaN(m,n);
  dtdz = NaN(m,n); 
  % Use multiple lookup tables to get times (handles variable shot depths) 
  % Cubic interpolation deals better with NaNs in lookup tables than spline
  for ilookup = 1:p.nlookup
    timePred1 = interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
              p.lookup(ilookup).ttTable,depth,range,'*linear',NaN);         
    timePred1 = reshape(timePred1,m,n);
    timePred(lookup==ilookup,:) = timePred1(lookup==ilookup,:);
    
    dtdr1 = interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
        p.lookup(ilookup).dtdrTable,depth,range,'*linear',NaN);   
    dtdr1 = reshape(dtdr1,m,n);
    dtdr(lookup==ilookup,:) = dtdr1(lookup==ilookup,:);
    
    dtdz1=interp2(p.lookup(ilookup).zTable,p.lookup(ilookup).rTable,...
        p.lookup(ilookup).dtdzTable,depth,range,'*linear',NaN);
    dtdz1 = reshape(dtdz1,m,n);
    dtdz(lookup==ilookup,:) = dtdz1(lookup==ilookup,:);
  end
  range = reshape(range,m,n);
  timePred(range>p.rangeMax) = NaN;

  %% Travel time misfits data vector
  d = time - timePred;
  d = d(:);
  d= d - tEvent(shotIndex,:);
  good = ~isnan(d);
  d = d(good);
  dataError = timeError(:);
  dataError = dataError(good);
  %range = range(good);
  %depth = depth(good);
  ntt = length(d);
  
  % Adjust data error
  rms(iter) = sqrt(mean(d.*d));
  rmsNorm(iter) = sqrt(mean(d.*d./dataError./dataError));
  if iter>1 && p.solve.adjustDataError
    dataError = dataError * rmsNorm(iter);
  end
    
  dMatrix = time - timePred;
  rmsEvent = NaN(m,1); rmsEventNorm = NaN(m,1);
  rmsStation = NaN(n,1); rmsStationNorm = NaN(n,1);
  for i = 1:m
    j = ~isnan(dMatrix(i,:));
    rmsEvent(i) = sqrt(mean(dMatrix(i,j).^2));
    rmsEventNorm(i) = sqrt(mean((dMatrix(i,j)./timeError(i,j)).^2));
  end
  for i = 1:n
    j = ~isnan(dMatrix(:,i));
    rmsStation(i) = sqrt(mean(dMatrix(j,i).^2));
    rmsStationNorm(i) = sqrt(mean((dMatrix(j,i)./timeError(j,i)).^2));
  end
      
  % Model misfits into data vector
  d = [d; p.srStation.x(ista)-xStation; p.srStation.y(ista)-yStation; zStation0-zStation; ...
          p.srEvent.x(ievt)-xEvent; p.srEvent.y(ievt)-yEvent; -tEvent;];
  % Apriori model uncertainties into data uncertainty vector
  dataError = [dataError; xStationErrorApr; yStationErrorApr; zStationErrorApr; ...
                          xEventErrorApr; yEventErrorApr; tEventErrorApr];
    
  %% Check for convergence 
  % Simple test based on minimal change in all model parameters and minimal
  % changes in rms and normalized rms over two iterations
  if iter>p.iterMax
    break
  end
  if iter>3 
    if all(abs(dm(~~ispace))<p.dxyzMax) && all(abs(dm(~ispace))<p.dtMax) && ...
       rms(iter-2)-rms(iter-1)<p.drmsMax && ...
       rms(iter-3)-rms(iter-2)<p.drmsMax  && ...
       rmsNorm(iter-2)-rmsNorm(iter-1) < p.drmsNormMax && ...
       rmsNorm(iter-3)-rmsNorm(iter-2) < p.drmsNormMax 
       break;
    end
  end  

  %% Calculate reciprocal of covariance matrix
  nmodel = 3*(n+m);
  Cr = sparse(1:nmodel+ntt,1:nmodel+ntt,1./dataError.^2,nmodel+ntt,nmodel+ntt);

  %% Spatial derivates for changes in locations
  azimuth = atan2(DX,DY);
  azimuth = reshape(azimuth,m,n);    
  dtdxStation = -dtdr.*sin(azimuth);
  dtdyStation = -dtdr.*cos(azimuth);
  dtdzStation = dtdz;
   
  %%Calculate A matrix (assume inverting for everything and then reduce)
  k=0; l=0;
  A=sparse([],[],[],ntt+nmodel,nmodel,6*n*m+nmodel);
  % Rows of A corresponding to data
  for i = 1:n
    for j = 1:m
      k = k+1;
      if good(k)
        l = l+1;    
        A(l,i) = dtdxStation(j,i);
        A(l,i+n) = dtdyStation(j,i);
        A(l,i+2*n) = dtdzStation(j,i);
        A(l,j+3*n) = -dtdxStation(j,i);
        A(l,j+3*n+m) = -dtdyStation(j,i);
        A(l,j+3*n+2*m) = 1;
      end
    end
  end
  %Rows of A corresponding to apriori model values
  for i = 1:n
    A(ntt+i,i) = 1;
    A(ntt+n+i,n+i) = 1;
    A(ntt+2*n+i,2*n+i) = 1;
  end
  for i = 1:m
    A(ntt+3*n+i,3*n+i) = 1;
    A(ntt+3*n+m+i,3*n+m+i) = 1;
    A(ntt+3*n+2*m+i,3*n+2*m+i) = 1;
  end
   
  % Columns of A that correspond to model parameters in inversion
  keep = logical(zeros(3*n+3*m+2,1));
  if p.solve.xyStation == 1
    keep(1:2*n) = 1;
  end
  if p.solve.zStation == 1
    keep(2*n+1:3*n) = 1;
  end
  if p.solve.xyEvent==1
    keep(3*n+1:3*n+2*m) = 1;
  end
  if p.solve.tEvent ==1
    keep(3*n+2*m+1:3*n+3*m) = 1;
  end

  % Select rows/columns of inversion components corresponding to parameters in inversion
  keep1 = logical([ones(ntt,1); keep]);
  A = A(:,keep);
  A = A(keep1,:);    
  d = d(keep1);    
  Cr = Cr(keep1,keep1); 

  % ispace is a logical which is true is a model parameter is spatial
  ispace = [ones(3*n+2*m,1); zeros(m,1)];
  ispace = ispace(keep);

  % Inversion (with aposteriori model Variances in Cm)
  if size(A,2)<p.nlsqr
    Cm = inv(A'*Cr*A);
    dm  = Cm *A'*Cr*d;  
    modelVar = full(diag(Cm));
  else
    [dm,flag] = lsqr(A'*Cr*A,A'*Cr*d,1e-6,1000);
    modelVar = zeros(size(dm));
    if flag
      disp('LSQR did not converge')
      %keyboard
    end
  end
   
  % Convert dm to changes in model parameters and modelVar to model errors
  k = 0;
  if p.solve.xyStation == 1
    dxStation = dm(1:n);
    dyStation = dm(n+1:2*n);
    xStationError = sqrt(modelVar(1:n));
    yStationError = sqrt(modelVar(n+1:2*n));
    k = k + 2*n;
  else
    dxStation = zeros(n,1);
    dyStation = zeros(n,1);
    xStationError = zeros(n,1);
    yStationError = zeros(n,1);
  end        
  if p.solve.zStation == 1
    dzStation = dm(k+1:k+n);
    zStationError = sqrt(modelVar(k+1:k+n));
    k = k + n;
  else
    dzStation = zeros(n,1);
    zStationError = zeros(n,1);
  end        
  if p.solve.xyEvent == 1
    dxEvent = dm(k+1:k+m);
    dyEvent = dm(k+m+1:k+2*m);
    xEventError = sqrt(modelVar(k+1:k+m));
    yEventError = sqrt(modelVar(k+m+1:k+2*m));
    k = k + 2*m;
  else
    dxEvent = zeros(m,1);
    dyEvent = zeros(m,1);
    xEventError = zeros(m,1);
    yEventError = zeros(m,1);
  end     
  if p.solve.tEvent == 1
    dtEvent = dm(k+1:k+m);        
    tEventError = sqrt(modelVar(k+1:k+m));        
    k = k + m;
  else
    dtEvent = zeros(m,1);
    tEventError = zeros(m,1);
  end  
  
  % Apply changes to Station and shot locations
  xStation = xStation + dxStation;
  yStation = yStation + dyStation;
  zStation = zStation + dzStation;
  xEvent = xEvent + dxEvent;
  yEvent = yEvent + dyEvent;
  tEvent = tEvent + dtEvent; 
  
  %% Update model Station depths if based on bathymetry
  if p.solve.useBathymetry
    [maplon,maplat] = xy2map(xStation,yStation,p.srGeometry);
    zStation0 = -interp2(p.srElevation.longitude,p.srElevation.latitude, ...
                     p.srElevation.data',maplon,maplat,'linear',NaN)/1000;
    if any(isnan(zStation0))
      disp('Problem getting Station depths form map')
      keyboard
    end
  end

end % End of inversion look


%% Outputs
s.converged = converged;
s.iter = iter;
s.xStation = xStation;
s.yStation = yStation;
s.zStation = zStation;
s.xStation_int = xStation_int;
s.yStation_int = yStation_int;
s.zStation_int = zStation_int;
s.xEvent_int = xEvent_int;
s.yEvent_int = yEvent_int;
s.tEvent_int = tEvent_int;
s.xEvent = xEvent;
s.yEvent = yEvent;
s.tEvent = tEvent;
s.xStationApr = p.srStation.x(ista);
s.yStationApr = p.srStation.y(ista);
s.zStationApr = zStation0;
s.xEventApr = p.srEvent.x(ievt);
s.yEventApr = p.srEvent.y(ievt);
s.tEventApr = zeros(size(tEvent));
s.xStationError = xStationError;
s.yStationError = yStationError;
s.zStationError = zStationError;
s.xEventError = xEventError;
s.yEventError = yEventError;
s.tEventError = tEventError;
s.rms = rms;
s.rmsNorm = rmsNorm;
s.rmsEvent = rmsEvent;
s.rmsEventNorm = rmsEventNorm;
s.rmsStation = rmsStation;
s.rmsStationNorm = rmsStationNorm;
s.ista = ista;
s.ievt = ievt;
s.time = time;
s.timePred = timePred;
s.timeError = timeError;

[mapx,mapy] = xy2map(xStation,yStation,p.srGeometry);
s.srStation.name = p.srStation.name(ista);
s.srStation.longitude = mapx;
s.srStation.latitude = mapy;
s.srStation.elevation = -zStation;
s.srStation.x = xStation;
s.srStation.y = yStation;
s.srStation.z = zStation;
if p.solve.xyStation
  s.srStation.xError = xStationError;
  s.srStation.yError = yStationError;
end
if p.solve.zStation
  s.srStation.elevationError = zStationError;
end
s.srStation.nsta = n;

[mapx,mapy] = xy2map(xEvent,yEvent,p.srGeometry);
s.srEvent.id = p.srEvent.id(ievt);
s.srEvent.latitude = mapy;
s.srEvent.longitude = mapx;
s.srEvent.type = p.srEvent.type(ievt);
s.srEvent.elevation = p.srEvent.elevation(ievt);
s.srEvent.x = xEvent;
s.srEvent.y = yEvent;
s.srEvent.z = p.srEvent.z(ievt);
s.srEvent.timeStatic = tEvent;
if p.solve.xyEvent 
  s.srEvent.xError = xEventError;
  s.srEvent.yError = yEventError;
end
if p.solve.tEvent
  s.srEvent.timeStaticError = tEventError;
end
s.srEvent.nevt = m;

%% plotting error ellipse  - Uncomment when plotting error ellipse
% 
% n1=[0;0;1];
% n2=[0;1;0];
% n3=[1;0;0];
% vp1=(eye(3)-n1*n1')*Cm*(eye(3)-n1*n1');
% vp2=(eye(3)-n2*n2')*Cm*(eye(3)-n2*n2');
% vp3=(eye(3)-n3*n3')*Cm*(eye(3)-n3*n3');
% 
% vp1=vp1(1:2,1:2);
% vp2=[vp2(1,1),vp2(1,3);vp2(3,1),vp2(3,3)];
% vp3=vp3(2:3,2:3);
% 
% [u1,lam1]=eig(vp1);
% [u2,lam2]=eig(vp2);
% [u3,lam3]=eig(vp3);
% 
% l1=sqrt(diag(lam1));
% l2=sqrt(diag(lam2));
% l3=sqrt(diag(lam3));
% dtheta=0.01;
% theta=0:dtheta:2*pi;
% 
% for i=1:length(theta)
%   v1=2.79*(l1(1)*u1(:,1)*cos(theta(i))+l1(2)*u1(:,2)*sin(theta(i)));
%   v2=2.79*(l2(1)*u2(:,1)*cos(theta(i))+l2(2)*u2(:,2)*sin(theta(i)));
%   v3=2.79*(l3(1)*u3(:,1)*cos(theta(i))+l3(2)*u3(:,2)*sin(theta(i)));
%   x1(i)=xStation(1)+v1(1);
%   y1(i)=yStation+v1(2);
%   x2(i)=xStation+v2(1);
%   z2(i)=zStation+v2(2);
%   y3(i)=yStation+v3(1);
%   z3(i)=zStation+v3(2);
%   
% end
% 
% s.a = x1;
% s.b = y1;

