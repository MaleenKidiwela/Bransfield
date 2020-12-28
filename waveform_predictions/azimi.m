function [t, pulse0, pulse, f, Qw, cw ] = azimi( N, Dt, x, c0, Q, f0 )
% input parameters:
% f0 corner frequency of Azimi Q model, in hz (e.g. 50)
% c0 base velocity in km/s (e.g. 4.5);
% x  propagation distance in km (e.g. 100)
% Q  low frequency quality factor (e.g. 10)
% N  number of samples in pulse (e.g. 1024);
% Dt sampling interval (e.g. 0.1)
% returned values
% t time array
% pulse0 input pulse, a unit spike at time N/2
% pulse attentated pulse
% f frequencies in Hz
% Qw frequency dependent quality factors
% cw frequency dependent phase velocities
% time series

t = Dt*[0:N-1]';
pulse0 = zeros(N,1);
pulse0(1)=1;
% standard fft setup
fny = 1/(2*Dt);
N2 = N/2+1;
df = fny / (N/2);
f = df*[0:N2-1]';
w = 2*pi*f;

w0 = 2*pi*f0;

% attenuation factor
% exp( -a(w) x ) = exp( - wx / 2Qc )
%
% propagation law with velocity c=w/k and slowness s=1/c=k/w % exp{ i(kx - wt) } = exp{ iw(sx - t) }
% propagation law
% exp( iwsx )
% Azimi's second law en.wikipedia.org/wiki/Azimi_Q_models %
% a(w) = a2 |w| / [ 1 + a3 |w| ] % note that for w<<w0 a(w) =
%
% s(w) = s0 + 2 a2 ln( a3 w ) / [ pi (1 - a3^2 w^2 ) ]
% now set a3 = 1/w0 where w0 is a reference frequency
% and set a2 = 1 / (2Qc0) where c0 is a reference velocity % so that
%a(w)=(1/2Qc0)|w|/[1+ |w/w0|]
% so for w/w0 << 1
% a(w) = w/(2Qc0) and Q(w) = w/(2 a c0)

a2 = 1 / (2*Q*c0);
a3 = 1 / w0;
a = a2*w ./ ( 1 + a3.*w );
Qw = w ./ (2.*a.*c0);
Qw(1) = Q;
ds = -2*a2*log(a3*w) ./ (pi*(1-(a3^2).*(w.^2 ))); ds(1)=0;
cw = 1./( (1/c0) + ds );
dt = fft(pulse0);
dp = dt(1:N2);
dp = dp .* exp(-a*x) .* exp(-complex(0,1)*w.*ds.*x);
dtnew = [dp(1:N2);conj(dp(N2-1:-1:2))]; % fold out negative frequencies
pulse = ifft(dtnew);

end

