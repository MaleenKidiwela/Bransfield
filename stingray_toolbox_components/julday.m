function [a,b]=julday(c,d,yr)
%Julian day conversion (both ways)
%If 1-2 input arguement
%  function [mo,day]=julday(jday,yr)
%If 2-3 input arguements
%  function [jday]=julday(mo,da,yr)
%  MO - Month (
%  DAY - Day
%  JDAY - Julian Day
%  YR - Year in '1995' format (must be greater than 1000) 
%       Default is a non-leap year
%  Vector inputs are okay

if nargin<1 | nargin>3 ;  error('JULDAY - Needs 1-3 input arguements'); end;
if nargin==1; forward=0; yr=1995; end
if nargin==2; if d(1)>1000; forward=0; yr=d; else forward=1; yr=1995; end; end
if nargin==3; forward=1; end; 
if nargin==2
  if length(c)==1 & length(d)>1 c=c*ones(size(d)); end;
  if length(d)==1 & length(c)>1 d=d*ones(size(c)); end;
  if length(c)~=length(d); error('JULDAY - Dimensions do not match'); end
end
if nargin==3
  if length(c)==1 & length(d)>1 c=c*ones(size(d)); end;
  if length(c)==1 & length(yr)>1 c=c*ones(size(yr)); end;
  if length(d)==1 & length(c)>1 d=d*ones(size(c)); end;
  if length(d)==1 & length(yr)>1 d=d*ones(size(yr)); end;
  if length(yr)==1 & length(c)>1 yr=yr*ones(size(c)); end;
  if length(yr)==1 & length(d)>1 yr=yr*ones(size(d)); end;
  if length(c)~=length(d) | length(c)~=length(yr); 
    error('JULDAY - Input dimensions do not match')
  end
end

for k=1:length(c);
  moday=[31 28 31 30 31 30 31 31 30 31 30 31];
  if (rem(yr(k),4)==0 & rem(yr(k),100)~=0) | rem(yr(k),400)==0
    moday(2)=moday(2)+1;
  end 
  moday=[0 cumsum(moday)];
  if forward
    a(k)=moday(c(k))+d(k);
  else
    a(k)=min(find(c(k)<=moday))-1;  
    b(k)=c(k)-moday(a(k));
  end
end
%end


