% one-off for display
%
%   j.a. collins

f = logspace(-1,2,2048);

% Hydrophone
[tf,tf2] = hightech_hti90u_whoi_obsip(f);   % volts/Pa
tf_hydrophone = tf(:);


amp = abs(tf_hydrophone);
faze = rad2deg(angle(tf_hydrophone));  

figure;
subplot(211);
loglog(f,amp); grid on;
axis([0.1 50 1e-5 3e-3]);
xlabel('Frequency (Hz)'); ylabel ('Amplitude (Volts/Pa)');

subplot(212);
semilogx(f,faze); grid on;
V = axis;
axis([0.1 50 -200 200]);
xlabel('Frequency (Hz)'); ylabel ('Phase (degrees)');

subplot(211)
title ('WHOI D2-OBS: High Tech HTI-90-U Hydrophone Response');


%%%% Q330
filter_name = 'FLinear-100.txt';
[H,b] = q330_firfilt_response(filter_name, f);

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


%%%% Combined Hydrophone and Quanterra Q330
tf = gain*tf_hydrophone.*H;

amp = abs(tf);
faze = rad2deg(angle(tf)); 

figure;
subplot(211);
loglog(f,amp); grid on;
axis([0.1 50 1 1e3]);
set (gca,'YTick',[1e0 1e1 1e2 1e3])
xlabel('Frequency (Hz)'); ylabel ('Amplitude (counts/Pa)');

subplot(212);
semilogx(f,faze); grid on;
V = axis;
axis([0.1 50 -180 180]);
xlabel('Frequency (Hz)'); ylabel ('Phase (degrees)');

subplot(211)
title ('WHOI D2-OBS: Combined High Tech HTI-90-U Hydrophone and Quanterra Q330 Response');