clear; close all

% Parameters
grx_db = 0;         % rx antenna gain
gtx_db = 0;         % tx antenna gain
ptx_dbw = 30;       % tx power, watts dB
freqc_hz = 3e9;     % carrier frequency
rcs_m2 = 1;         % radar cross section, meter sq
bw_hz = 1e6;        % rx noise bandwidth
noisefig_db = 10;   % system noise figure

range_mi = 1:100;   % one way range to target, miles
temp = 290;         % standard temperature, kelvin
k = 1.38e-23;       % Boltzman's constant
c = 299702458;      % speed of light

range_m = 1609.34*range_mi;

snrrx_dbm = ptx_dbw + gtx_db + grx_db + 10*log10(rcs_m2) + ...
    20*log10(c/freqc_hz) - 10*3*log10(4*pi) - 10*4*log10(range_m) - ...
    10*log10(k*temp) - 10*log10(bw_hz) - noisefig_db;

plot(range_mi, snrrx_dbm)
ylabel('SNR (dBW)')
xlabel('One-way Range (miles)')
grid; grid minor

