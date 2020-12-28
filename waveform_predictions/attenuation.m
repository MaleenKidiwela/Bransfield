clear all;

%% Attenuation from ground

Q1 = 100;  % quality factor at low frequencies
f0 = 200;  % frequency below which
N=1000;    % number of samples
Dt=0.005; % sampling interval
c0=5;      % low frequency velocity

x=100;  % propagation distance
[ t, pulse0, pulse1, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=50;
[ t, pulse0, pulse2, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=30;
[ t, pulse0, pulse3, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=25;
[ t, pulse0, pulse4, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=20;
[ t, pulse0, pulse5, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=15;
[ t, pulse0, pulse6, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=10;
[ t, pulse0, pulse7, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

x=5;
[ t, pulse0, pulse8, f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );

% plot Azimi pulses
figure(1);
clf;
hold on;
plot( t, pulse0, 'k-', 'LineWidth', 2 );
plot( t, pulse1,'LineWidth', 1 );
plot( t, pulse2,'LineWidth', 1 );
plot( t, pulse3,'LineWidth', 1 );
plot( t, pulse4,'LineWidth', 1 );
plot( t, pulse5,'LineWidth', 1 );
plot( t, pulse6,'LineWidth', 1 );
plot( t, pulse7,'LineWidth', 1 );
plot( t, pulse8,'LineWidth', 1 );
legend('100','50','30','25','20','15','10','5')
title('azimi pulse for Q=100'); 
xlabel('t');
ylabel('u');


%% Q330 Filter response

dt = 1/2000;
imp = zeros(1000,1);
imp(1)=1;

kpts = length(imp);
dff = 1/(kpts*dt);
ff_nyq = 1/(2*dt);
ff = [0:dff:ff_nyq];
ff = ff(:);

% filter impulse with Q330 minimum-phase filter
[H,b,delay] = q330_firfilt_response('FLinear-200.txt', ff);
st1 = filter(b,1,imp);

%% Geophone response

dt = 1/2000;
imp = zeros(1000,1);
imp(1)=1;

kpts = length(imp);
dff = 1/(kpts*dt);
ff_nyq = 1/(2*dt);
ff = [0:dff:ff_nyq];
ff = ff(:);

[tf,tf2] = geospace_gs11d_whoi_obsip(imp,'v');   % volts/m/s
tf_geophone = tf(:);


amp = abs(tf_geophone);

%% airgun

load('airgun.mat')

%% convolve

conv1= conv(airgun(:,2),pulse3);
conv2 = conv(conv1,st1);
conv3 = conv(conv2,amp);

plot(real(conv3))