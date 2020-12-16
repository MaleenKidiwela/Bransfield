function [H,b,delay] = q330_firfilt_response(filter_name, f);
%
%   Function "q330_firfilt_response.m" generates the composite frequency 
% response of the Quanterra Q330 FIR filters.  Output is dimensionless.
% The variable "f" is a frequency vector.    The variable "filter_name" 
% is one of:
%
% FLbelow100-1.txt   FLbelow20-1.txt    FLbelow40-1.txt    FLinear-1.txt
% FLbelow100-10.txt  FLbelow20-10.txt   FLbelow40-10.txt   FLinear-10.txt
% FLbelow100-100.txt FLbelow20-100.txt  FLbelow40-100.txt  FLinear-100.txt
% FLbelow100-20.txt  FLbelow20-20.txt   FLbelow40-20.txt   FLinear-20.txt
% FLbelow100-200.txt FLbelow20-200.txt  FLbelow40-200.txt  FLinear-200.txt
% FLbelow100-40.txt  FLbelow20-40.txt   FLbelow40-40.txt   FLinear-40.txt
% FLbelow100-50.txt  FLbelow20-50.txt   FLbelow40-50.txt   FLinear-50.txt
%
% See documents:
% "Q330 Response-v1_0.pdf" by Joe Steim;
% Also, within MATLAB, type "help filter".
%
% To see poles and zeros: 
% a=1; [z,p,k] = tf2zpk(b,a); zplane(z,p);
%
% In the old days (pre July 2007), coefficients in the SEED response file were 
% reversed relative to the order of the coefficients in Steim's file.  
% Appendix C of SEED Manual says: 
% "In July 2007, the FDSN adopted a convention that requires the coefficients to be 
% listed in forward order as used in equation 8.  As a reference, minimum-phase filters 
% (which are asymmetric) should be written with the largest values near the beginning 
% of the coefficient list.
%
% Here, we use the unreversed (time-forward) coefficients.
% See report by Bob Herrmann that shows that "evalresp" (V3.2.37) generates 
% an incorrect response. 
%
% Note that to check the phase delay of this filter, either use the "b" 
% vector or the uncorrected "H" (i.e without remiving the delay) response.  Then:
% faze = unwrap(angle(H)); faze_delay = -faze./(2*pi*f); semilogx(f,faze_delay); grid on;  
%
% To check the impulse response do:
% imp = [1; zeros(256,1)]; a = 1;
% y = filter(b,a,imp); plot(y(1:50)); grid on;
% OR
% pltseis(y,dt,-delay,0,1);  % dt is appropriate for filter selected
% OR
% fvtool(b,a);
%
% USAGE: [H,b,delay] = q330_firfilt_response(filter_name, f);
%                                                            j.a.collins
%-----------------------------------------------------------------------

disp(' ');

fname_base = '/Users/jac/WHOI_OBS/Quanterra_Stuff/Q330/Documents/Q330Response/Q330FIR-Response';
fname_base = './';
fname = strcat(fname_base,'/',filter_name);
fprintf(1, '    Filter file: %s\n', fname);
fid = fopen(fname,'r');
if (fid == -1)
    error (['    Cannot open file: ', fname]);
end


jnk = textscan(fid,'%s %*[^\n\r]',1);
filter_name = textscan(fid,'%s%*[^\n]',1);
filter_delay = textscan(fid,'%f%*[^\n]',1);
filter_length = textscan(fid,'%f%*[^\n]',1);
filter_gain = textscan(fid,'%f%*[^\n]',1);
filter_output_rate = textscan(fid,'%f%*[^\n]',1);
jnk = textscan(fid,'%s %*[^\n\r]',1);
jnk = textscan(fid,'%s %*[^\n\r]',1);
filter = textscan(fid,'%f',filter_length{1});
fclose(fid);

fprintf(1, '    Filter Name: %s\n',char(filter_name{1}));
fprintf(1, '    Filter Delay (s): %f\n',filter_delay{1});
fprintf(1, '    Filter Length (points): %d\n',filter_length{1});
fprintf(1, '    Filter Gain: %f\n',filter_gain{1});
fprintf(1, '    Filter Output Rate: %d\n',filter_output_rate{1});

Fs = filter_output_rate{1};
composite_delay = filter_delay{1};
npts = filter_length{1};
b = filter{1};
if (npts ~= length(b))
    disp ('    Error reading filter coefficients!');
end

f = f(:);
b = b(:);

dt = 1/Fs;
f_nyq = Fs/2.0;
ndx = f <= f_nyq;
%%f = f(ndx);  % use this to restrict output to frequencies <= Nyquist

H = zeros(size(f));
for n = 1:length(f)
    tmp = 0.0;
    for nn = 1:length(b);
        H(n) = b(nn)*exp(-i*2*pi*(nn-1)*f(n)*dt) + tmp;
        tmp = H(n);
    end
end

% The time shift "composite_delay" introduced at the Q330 is 
% removed by the Data Processor (e.g., Baler or Pecos or Antelope).
% Note: Antelope 4.8 and DMC use the "q330_delay".
%%q330_delay = composite_delay - dt;
%%delay = q330_delay;

% subtract delay;  comment out the following 2 lines if checking the pahse delay
delay = composite_delay;       % this is correct
H = H.*exp(i*2*pi*f*delay);     



H = H(:);


return
