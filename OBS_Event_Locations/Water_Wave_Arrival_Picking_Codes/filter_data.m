%----Butterworth filter
if ~length(fmini); fmini=0; end;
if ~length(fmax); fmax=0; end;
if fmini
  wmin=max(0,min(1,abs(fmini*sampint*2)));
end
if fmax
  wmax=max(0,min(.99,abs(fmax*sampint*2)));
end
if fmini & ~fmax
  db=dboct;
  if dboct==0; db=24; end
  o=buttord(wmin,wmin/2,3,abs(db));
  if db>0
    [b,a]=butter(o,wmin,'high');
    seis=filter(b,a,seis_unfilt);
  else
    [b,a]=butter(ceil(o/2),wmin,'high');
    seis=filtfilt(b,a,seis_unfilt);
  end
  disp(['  Filtering data with a low-cut filter at ' num2str(fmini) ' Hz ']); drawnow
elseif fmax & ~fmini
  db=dboct;
  if dboct==0; db=24; end
  o=buttord(wmax,wmax/2,3,abs(db));
  if db>0
    [b,a]=butter(o,wmax);
    seis=filter(b,a,seis_unfilt);
  else
    [b,a]=butter(ceil(o/2),wmax);
    seis=filtfilt(b,a,seis_unfilt);
  end
  disp(['  Filtering data with a high-cut filter at ' num2str(fmax) ' Hz ']);
elseif fmax & fmini
  db=dboct;
  if dboct==0; db=24; end
  o=buttord(wmin,wmin/2,3,abs(db));
  if db>0
    [b,a]=butter(o,[wmin wmax]);
    seis=filter(b,a,seis_unfilt);
  else
    [b,a]=butter(ceil(o/2),[wmin wmax]);
    seis=filtfilt(b,a,seis_unfilt);
  end
  disp(['  Filtering data with a bandpass filter at ' num2str(fmini) '-' num2str(fmax) ' Hz ']); drawnow
elseif ~fmax & ~fmini
  seis=seis_unfilt;
  disp(['  Data is now unfiltered ']); drawnow
end
% Remove filter buffer
seis=seis(nfilterbuf:end-nfilterbuf+1,:);

% MUTE NO LONGER NEEDED WITH FILTER BUFFER OPTION
%----Mute start (+ end for zero phase filters) of SEIS columns 
%if fmini; fmute=fmini; else; fmute=fmax; end;
%if fmute
%  nmute=ceil((1/fmute)/sampint(1));
%  seis(1:nmute,:)=zeros(nmute,n);
%  ramp=linspace(0,1,nmute)';
%  ramp=ramp(:,ones(n,1));
%  seis(nmute+1:2*nmute,:)=seis(nmute+1:2*nmute,:).*ramp;
%  if db<0
%    seis(end-nmute+1:end,:)=zeros(nmute,n);
%    ramp=flipud(ramp);
%    seis(end-2*nmute+1:end-nmute,:) = ...
%               seis(end-2*nmute+1:end-nmute,:).*ramp;
%  end
%  disp('  SEGY data successfully filtered')
%  disp(' '); drawnow
%end
