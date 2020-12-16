function [tf,tf2] = geospace_gs11d_whoi_obsip(f,units);
%
%   Function "geospace_gs11d_whoi_obsip" calculates the transfer function between ground
% velocity (m/s) and digital counts as a function of frequency (f) for the
% Geospace GS-11D geophone.   
%
% Set variable "units" = 'v' for velocity output (default)
% or "units" = 'a' for acceleration.

%   The vector "tf" holds the complex-valued "frequency response function", and the vector
% "tf2" hold the square of the absolute value of the "frequency response function".
%   References: (i) Scherbaum, Chapter 4.
%
% Usage:  [tf,tf2] = geospace_gs11d_whoi_obsip(f,units);
%                                                         ---j.a.collins, whoi 
% ****************************************************************************

if ((nargin < 2) | isempty(units))
    units = 'v';
end
cm_per_inch = 2.54;

clear p z
nzeros = 2;
npoles = 2;



%########## Stage 1 #############

Mass = 23.6e-3;                         % from data sheet; kilograms
damping_constant_mechanical = 0.35;     % from data sheet

EDC = 2.54;                             % Electro-dynamic constant (V/inch/s)
EDC = (EDC/cm_per_inch) * 1e02;         % V/m/s

Rc = 4000;  % Coil resistance (ohms);
Rs = 18200; % Shunt resistance (ohms);
sensitivity = EDC*Rs/(Rs+Rc);


f0 = 4.5;                               % Resonant frequency of seismometer (Hz)
w0 = 2*pi*f0;                           % 

h = 0.7;                                % usual asssumption
damping_constant_electrical = EDC^2/(2*Mass*w0*(Rc + Rs));
h = damping_constant_mechanical + damping_constant_electrical;


p(1) = -( h + sqrt(h^2 - 1) )*w0;
p(2) = -( h - sqrt(h^2 - 1) )*w0;

z(1) = complex(0.0,0.0);
z(2) = complex(0.0,0.0);

f = f(:);
w = 2*pi*f;
s = i*w;
Hn = ones(size(s));
Hd = ones(size(s));
for n = 1:nzeros
    Hn = Hn.*(s-z(n));
end
for n = 1:npoles
    Hd = Hd.*(s-p(n));
end
H = Hn./Hd;


% find value of transfer function at a normalization frequency of 20 Hz
f_norm = 20;
s_norm = i*2*pi*f_norm;
Hn_norm = 1;
Hd_norm = 1;
for n = 1:nzeros
    Hn_norm = Hn_norm*(s_norm-z(n));
end
for n = 1:npoles
    Hd_norm = Hd_norm*(s_norm-p(n));
end
H_norm = abs(Hn_norm/Hd_norm);


normalization = 1/H_norm;
A0 = normalization;
G = sensitivity/A0;
tf = G*A0*H;  % units are volts/m/s
tf = tf(:);
tf2 = abs(tf).^2;                                 % units are V**2/(m/s)**2
%

disp(' ');
fprintf(1, 'Sensitivity (V/(m/s)): %10.4E\n', G);
disp(' ');

disp(' ');
fprintf(1, 'Normalization: %10.4E\n', A0);
fprintf(1, 'Normalization Frequency (Hz): %10.4E\n', f_norm);
disp(' ');

disp(' ');

fprintf(1, 'Number of Poles: %d\n', npoles);
fprintf(1, 'Poles (radians/second): \n');
for n = 1:npoles
    fprintf (1, '%s%10.4E, %10.4E%s\n','complex(',real(p(n)),imag(p(n)),')');
end    
disp(' ');

fprintf(1, 'Number of Zeros: %d\n', nzeros);
fprintf(1, 'Zeros (radians/second): \n');
for n = 1:nzeros
    fprintf (1, '%s%10.4E, %10.4E%s\n','complex(',real(z(n)),imag(z(n)),')');
end    
disp(' ');

if (strncmp(lower(units),'v',1)) 
    return;
end
if (strncmp(lower(units),'a',1)) 
    tf = tf./(s + eps);      % units are V/m/s^2
end
tf = tf(:);
tf2 = abs(tf).^2;              % units are V**2/m**2/s**4

return;