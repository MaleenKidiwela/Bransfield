function [secnds] = date2secnds(year,day,hour,minute,second);
%
% total # of seconds since 00:00:00 (UTC) Jan 1, 1970 
% Function date2secnds finds total # of seconds since 
% 00:00:00 (UTC) Jan 1, 1970, given year, day of year, hour, 
% minute, and seconds (includes leap seconds)
%
% USAGE: [secnds] = date2secnds(year,day,hour,minute,second);
%                                                            j.a.collins
%-----------------------------------------------------------------------
%
% Leap year bug fixed by William Wilcock 11/8/96

yr_ref = 1970;
day_ref = 1; 
yr2sec = 31536000;
dy2sec = 86400;
hr2sec = 3600;
mn2sec = 60;

secnds = 0.0;

delta_yrs = year - yr_ref;
for n = 1:delta_yrs
  secnds = secnds + yr2sec;
  yr = yr_ref + n - 1;
% add leap year 
  leap = rem(yr,4) == 0 & rem(yr,100) ~= 0 | rem(yr,400) == 0; 
    if (leap == 1) 
      secnds = secnds + dy2sec;
    end
end
secnds = secnds + (day-1)*dy2sec + hour*hr2sec + minute*mn2sec + second;
%
% add leap seconds
%[nleap] = leap_secnds (secs);
%secnds = secnds + nleap;
        
