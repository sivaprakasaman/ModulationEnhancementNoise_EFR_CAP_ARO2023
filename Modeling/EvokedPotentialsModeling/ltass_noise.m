function noise = ltass_noise(fs,SPL,varargin)
%ltass_noise: Generate a vector of LTASS noise.
% Syntax:
%   NOISE = ltass_noise(Fs,SPL,n,m)
% where Fs is the sample rate in samples/sec, SPL is the desired sound
% pressure level of the noise in decibels (dB SPL), n and m are the number
% of rows and columns respectively of the noise matrix.  NOISE is a matrix
% of Long Term Average Speech Spectrum (LTASS) noise with the desired level
% and size.

% Long Term Average Speech Spectrum (LTASS) noise is random noise that has
% the same spectrum as a long-term average of real speech.  It is useful
% as a masking noise in speech intelligibility experiments.

% Written by Douglas M. Schwarz
% douglas.schwarz@rochester.edu
% 16 December 2014


% ---------------------------------------------------------------------
% Reference:
%
% Byrne, Denis, et al. "An international comparison of long-term average
% speech spectra." The Journal of the Acoustical Society of America 96.4
% (1994): 2108-2120.
% ---------------------------------------------------------------------


fn = fs/2; % Nyquist rate
Pref = 20e-6; % reference pressure in pascals (Pa)
desired_rms = Pref*10.^(SPL/20);

% LTASS noise magnitude specification.
M_dB = [30 38.6 54.4 57.7 56.8 60.2 60.3 59.0 62.1 62.1 60.5 56.8 53.7 ...
	53.0 52.0 48.7 48.1 46.8 45.6 44.5 44.3 43.7 43.4  41.3  40.7 30].';
f = [0 80   100  125  160  200  250  315  400 500  630  800  1000 ...
	1250 1600 2000 2500 3150 4000 5000 6300 8000 10000 12500 16000 25e3].';

% Build a filter with freq. response equal to the LTASS specification.
use = f > 0 & f < fn;
f2 = [0;f(use);fn];
M_dB2 = interp1(f,M_dB,f2,'linear','extrap');
M2 = 10.^(M_dB2/20);
b0 = fir2(5000,f2/fn,M2);

% % Apply filter to white noise.
% noise = conv(randn(varargin{:}),b0,'same');
% 
% % Scale noise to desired SPL.
% noise = noise*(desired_rms/rms(noise));


% Apply filter to white noise.
noise = conv(randn(varargin{:}),b0,'same');

% Scale noise to desired SPL (integral method).
[H,ff] = freqz(b0,1,fs,fs);
rms_infinity = sqrt(2*trapz(ff,abs(H).^2)/fs);
noise = noise*(desired_rms/rms_infinity);
