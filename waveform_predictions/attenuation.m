clear all;

%% Attenuation from ground

Q1 = 100;   % quality factor at low frequencies
f0 = 2500;  % frequency below which
N = 10000;  % number of samples in azimi's time series
N1 = 1000;  % Number of samples in geophone and quanterra responses 
Dt = 0.0005; % sampling interva
c0=5;      % low frequency velocity

i=0;
for x=[30,25,20,15,10,5,1]
    i=i+1;
    [ t, pulse0, pulse(:,i), f, Qf1, cw1] = azimi( N, Dt, x, c0, Q1, f0 );
    pulse(:,i) = real(pulse(:,i));
end


%% Q330 Filter response

imp = zeros(N1,1);
imp(1)=1;

kpts = length(imp);
dff = 1/(kpts*Dt);
ff_nyq = 1/(2*Dt);
ff = [0:dff:ff_nyq];
ff = ff(:);

% filter impulse with Q330 minimum-phase filter
[H,b,delay] = q330_firfilt_response('FLinear-200.txt', ff);
st1 = filter(b,1,imp);

figure(1)
subplot(3,1,1)
plot(st1)
title('Q330 Response'); 
hold on
%% Geophone response

% dt and imp should only be defined in one place so deleted from here

kpts = length(imp);
dff = 1/(kpts*Dt);
ff_nyq = 1/(2*Dt);
ff = [0:dff:ff_nyq];
ff = ff(:);

tf = geospace_gs11d_whoi_obsip(ff,'v');   % volts/m/s
tf(end) = real(tf(end));
tf = [tf; conj(tf(end-1:-1:2))];
st2 = real(ifft(tf));

figure(1)
subplot(3,1,2)
plot(st2)
title('Geophone Response'); 

%% airgun

load('airgun.mat')

%% convolve

for i= 1:7
    conv1= conv(airgun(:,2),pulse(:,i));
    conv2 = conv(conv1,st1);
    conv3(:,i) = conv(conv2,st2);
end

figure(2)
clf
hold on
for i= 1:7
    plot((0:length(conv3(:,i))-1)*Dt,real(conv3(:,i)))
end
xlim([0 1])
legend('x=30','x=25','x=20','x=15','x=10','x=5','x=1')
title('Waveforms predicted')
xlabel('t')

%%
% plot Azimi pulses
figure(1)
subplot(3,1,3)
hold on
plot( t, pulse0, 'k-', 'LineWidth', 2 );
for i=1:7
    plot( t, pulse(:,i),'LineWidth', 1 );
end

axis([0,0.5,0,0.16])
legend('pulse0','x=30','x=25','x=20','x=15','x=10','x=5','x=1')
title('Azimi pulse for Q=100 varying ranges'); 
xlabel('t');
ylabel('u');


%% applying Filter

load('currentMenu.mat')
currentMenu=currentMenuFilt;
traceMetaData.samprate(1)=2000;
constantSamprate=true;
flim(1) = currentMenu.filter.lim0;
flim(2) = currentMenu.filter.lim1;
wlim = flim * 2 * (1/traceMetaData.samprate(1));


for i=1:7
    traceDataUnfilt=conv3(:,i);
    %butter_filterSEGY
    [b,a] = butter(ceil(currentMenu.filter.order),wlim(1),'high');
    traceData = filter(b,a,traceDataUnfilt);
    
    filtered_data(:,i)= traceData;
    clear traceData b a traceDataUnfilt
    disp(['Filtered with a band-pass filter at [' num2str(flim(1)) '] Hz '])
end

figure(3)
clf
hold on
for i= 1:7
    plot((0:length(filtered_data(:,i))-1)*Dt,real(filtered_data(:,i)))
end

xlim([0 1])
legend('x=30','x=25','x=20','x=15','x=10','x=5','x=1')
title('Waveforms predicted filterd')
xlabel('t')

