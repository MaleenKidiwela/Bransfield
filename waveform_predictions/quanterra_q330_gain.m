function [gain] = quanterra_q330_gain();
%
%   Function "q330_gain.m" returns the gain of the Quanterra Q330 
% "analog-to-digital" converter. The gain represents the digital 
% counts output by the the digitizer as a function of input volts.
%
%
% USAGE: [gain] = quanterra_q330_gain;
%                                                            j.a.collins
%-----------------------------------------------------------------------

gain = 419430;  % (counts/V)
return;