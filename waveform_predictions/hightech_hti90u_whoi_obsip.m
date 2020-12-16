function [tf,tf2] = hightech_hti90u_whoi_obsip(f);
%
%   Function "hightech_hti90u_whoi_obsip" calculates the transfer function (V/Pa) between 
% pressure and volts as a function of frequency (f) for the High Tech Inc, HTI-90-U 
% hydrophone with built-in preaamplifier.
%
%   References: (i) "3HightechTransferFn.pdf" by K. Peal, 23 May 2007; 
% (ii) Scherbaum, Chapter 4.
%
% Appendix C of SEED Manual, eqn 6:  At any frequency f (in Hz) the response is:
%    G(f) = Sd*A0*(Prod(s - zn)/Prod(s-pn)) = Sd*A0*Hp(s)
%
% where s = i*2*pi*f if the reference frequency is 1 radian/second, 
% and s =i*f if the reference frequency is 1 Hz.
%
% Appendix C of SEED manual says to partition the response by choosing A0 
% so that the modulus of A0 times the modulus of the ratio of polynomials 
% equals 1.0 at the normalizing frequency fn; the Sd specified in Blockette [58] 
% is then the stage gain at that frequency, so |G(fn)| = Sd.
%
% ####################################################################################
%  Added March 01, 2017.  Emily Hooft pointed out that the Santorini data showed that
% the hydrophone shows a clear negative signal when the Langseth shoots from near-zero
% offset. Tim put one of our hydrophones in his DPG test gear and confirmed that this
% is so.  
%
% From Brian S. Spychalski <hightechinc@att.net>
% "The geophysical convention is a pressure decrease gives a positive output and a pressure 
% increase gives a negative output.  So it sounds like your phones are operating correctly." 
% ####################################################################################
%  
%
% Usage:  [tf,tf2] = hightech_hti90u_whoi_obsip(f);
%                                                         ---j.a.collins, whoi 
% ****************************************************************************

w = 2*pi*f;
s = i*w;
clear p z


C_1   = 0.1e-6;    % Farads
C_H   = 800e-12;   % Farads;
C_add = 47e-6;     % Farads;

R_1   = 100e3;     % Ohms
R_2   = 100e3;     % Ohms
R_3   = 49.9e6;    % Ohms
R_6   = 4.99e3;    % Ohms
R_8   = 10e3;      % Ohms
R_add = 10e3;      % Ohms

%########## Hyrophone Element ##########
sensitivity_hydrophone = 5.0119e-10;     % V/micro_Pa or -186 dB re. 1V/micro_Pa
sensitivity_hydrophone = sensitivity_hydrophone*1e6;    % V/Pa

%########## Hydrophone Preamp #############
R_E2 = R_1*R_2/(R_1 + R_2);
C_E1 = C_1*C_H/(C_1 + C_H);
a_3 = R_3*R_E2*C_E1*C_add;
b_3 = R_E2*C_add + R_E2*C_E1 + R_3*C_E1 + (R_3*R_E2*C_E1)/R_6; 
c_3 = 1;

num = R_3*(R_8 + R_add)*C_E1.*s.*(R_E2*C_add*s + 1);
denom = R_add*(a_3*s.^2 + b_3*s + c_3);
% =>
% 2 zeros at: (1) s = 0; (2) s = -1/(R_E2*C_add)
% 2 poles at: (-b_3 +/- sqrt(b_3^2 - 4*a_3*1))/(2*a_3)

num_coef = R_3*(R_8 + R_add)*C_E1*[R_E2*C_add 1 0];
denom_coef = R_add*[a_3 b_3 c_3];
[z,p,k] = tf2zp(num_coef,denom_coef);

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% Here we account for the fact that a pressure decrease gives a positive output.
% This -ve sign added on March 01, 2017.  Calculations of the transfer function
% prior to this date are in error!!!!!!!!!!!!!!!
k = -1*k; 
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nzeros = length(z);
npoles = length(p);

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
H = Hn./Hd;   % yes, do not include gain 'k'.


% find value of transfer function at a normalization frequency of 500 Hz
f_norm = 500;
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
G = sensitivity_hydrophone*k/A0;    % remember to include gain 'k'.
tf = G*A0*H;                        % units are V/Pa
tf2 = abs(tf).^2;                   % units are V**2/Pa**2
%


disp(' ');
fprintf(1, 'Sensitivity (V/Pa): %10.4E\n', G);
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

keyboard
return;