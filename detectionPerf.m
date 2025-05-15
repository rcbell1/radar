clear; close all

% Detection performance curve comparisons
n = 1000; % number of samples
snr_db = -30:0.2:30;
pfa = 1e-5;

snr = 10.^(snr_db/10);

%% Matched filter
pd_mf = qfunc(qfuncinv(pfa) - sqrt(n*snr));

%% Energy Detector Exact
thresholds = chi2inv(1 - pfa, 2*n);
pd_ed = 1 - ncx2cdf(thresholds, 2*n, 2*n*snr);

%% N>>1 Energy Detector
pd_ed1 = qfunc(qfuncinv(pfa) - sqrt(n)*snr);

%% N = 1 Radar Energy Detector
pd_ed2 = pfa.^(1./(1+snr));

%% Plots
semilogy(snr_db, pd_mf, 'k-'); hold all
semilogy(snr_db, pd_ed, '-.')
semilogy(snr_db, pd_ed1, ':')
semilogy(snr_db, pd_ed2, '--')
grid on; grid minor
legend('MF', 'Exact ED', 'N>>1 ED', 'Radar Approx N=1 ED', 'Location', 'northwest')

