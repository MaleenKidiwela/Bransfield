%  filter impulse with Q330 minimum-phase filter to evaluate time of first-break.

% make impulse
dt = 1/200;
imp = [zeros(500,1);1;zeros(523,1);zeros(1014,1)];

kpts = length(imp);
dff = 1/(kpts*dt);
ff_nyq = 1/(2*dt);
ff = [0:dff:ff_nyq];
ff = ff(:);

% filter impulse with Q330 minimum-phase filter
[H,b,delay] = q330_firfilt_response('FLinear-200.txt', ff);
st = filter(b,1,imp);

% subtract delay
%delay = round(delay/dt)*dt;  % this avoids some distortion.
stff = fft(st);
tmp0 = stff(1:kpts/2 + 1) .* exp(i*2*pi*ff*delay);
tmp1 = [tmp0; conj(flipud(tmp0(2:end-1)))];
stf = real(ifft(tmp1));


% plot
figure
subplot(311);
pltseis(imp,dt,0,2.4,2.7,'bo-'); grid on
subplot(312);
pltseis(st,dt,0,2.4,2.7,'bo-'); grid on;
subplot(313);
pltseis(stf,dt,0,2.4,2.7,'bo-'); grid on;

figure
subplot(311);
pltseis(imp,dt,0,2.4,2.7,'bo-'); grid on
subplot(312);
pltseis(st,dt,0,2.4,2.7,'bo-'); grid on;
subplot(313);
pltseis(stf,dt,0,2.4,2.7,'bo-'); grid on;
