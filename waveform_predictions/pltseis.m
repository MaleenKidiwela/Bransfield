function [ylim] = pltseis(seis,dt,t0,tst,tnd,line,scale);
%
% Usage: [ylim] = pltseis(seis,dt,t0,tst,tnd,line,scale);
% e.g. [ylim] = pltseis(seis,dt,t0,tst,tnd,'m-','n');

set_hold = 0;
if (ishold)
    set_hold = 1;
end

[mpts,nseis] = size(seis);
tt = t0 + [0:dt:(mpts-1)*dt];

if (nargin < 7)
  scale = 't';
end
if (nargin < 6)
  line = '-';
end

if (tst == 0 & tnd == 0)
    tst = 0;
    tnd = (mpts-1)*dt;
end
nst = max(1,round((tst-t0)/dt) + 1);
nst = min(nst,mpts);
nnd = min(mpts,round((tnd-t0)/dt));
nnd = max(nnd,1);


ylim(1) = min(min(seis(nst:nnd,1:nseis)));
ylim(2) = max(max(seis(nst:nnd,1:nseis)));
if (scale == 'n')
  smax = max(max(abs(seis(nst:nnd,1:nseis))));
else
  smax = 1;
end
ylim(1) = ylim(1)/smax;
ylim(2) = ylim(2)/smax;
for n=1:nseis
  plot (tt(nst:nnd),seis(nst:nnd,n)/smax,line)
  hold on
end
xlabel ('Time (s)')
ylabel ('Amplitude')
grid on;

axlim = axis;
axis([tst tnd axlim(3) axlim(4)]);

if (set_hold)
    hold on;
else
    hold off;
end

return
