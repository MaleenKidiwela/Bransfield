% one-off for display
%
%   j.a. collins


f = logspace(-1,2,2048);

%%%% Geophone
[tf,tf2] = geospace_gs11d_whoi_obsip(f,'v');   % volts/m/s
tf_geophone = tf(:);


amp = abs(tf_geophone);
faze = rad2deg(angle(tf_geophone));  

figure;
subplot(211);
loglog(f,amp); grid on;
axis([0.1 50 1e-2 3e2]);
set (gca,'YTick',[1e-2 1e-1 1e0 1e1 1e2])
xlabel('Frequency (Hz)'); ylabel ('Amplitude (Volts/m/s)');

subplot(212);
semilogx(f,faze); grid on;
V = axis;
axis([0.1 50 -200 200]);
xlabel('Frequency (Hz)'); ylabel ('Phase (degrees)');

subplot(211)
title ('WHOI D2-OBS: Geospace GS-11D Geophone Response');


%%%% Q330
filter_name = 'FLinear-100.txt';
[H,b] = q330_firfilt_response(filter_name, f);

preamp_gain = 30;
[gain] = quanterra_q330_gain;   % counts/volt
tf = gain.*H;

amp = abs(tf);
faze = rad2deg(angle(tf));  

figure;
subplot(211);
loglog(f,amp); grid on;
axis([0.1 50 1e5 1e6]);
xlabel('Frequency (Hz)'); ylabel ('Amplitude (counts/volt)');

subplot(212);
semilogx(f,faze); grid on;
V = axis;
axis([0.1 50 -45 45]);
xlabel('Frequency (Hz)'); ylabel ('Phase (degrees)');

subplot(211)
title ('WHOI D2-OBS: Quanterra Q330 FLinear-100 Response');




%%%% Combined Geophone and Quanterra Q330
tf = gain*preamp_gain*tf_geophone.*H;

amp = abs(tf);
faze = rad2deg(angle(tf)); 

figure;
subplot(211);
loglog(f,amp); grid on;
axis([0.1 50 1e5 1e10]);
set (gca,'YTick',[1e5 1e6 1e7 1e8 1e9 1e10])
xlabel('Frequency (Hz)'); ylabel ('Amplitude (counts/m/s)');

subplot(212);
semilogx(f,faze); grid on;
V = axis;
axis([0.1 50 -200 200]);
xlabel('Frequency (Hz)'); ylabel ('Phase (degrees)');

subplot(211)
title ('WHOI D2-OBS: Combined Geospace GS-11D Geophone and Quanterra Q330 Response');