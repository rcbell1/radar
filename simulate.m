close all; clear
rng(10024)

%% Transceiver parameters
sampleRate = 100e6;         % Sample rate (Hz)
prf = 1e3;                  % Pulse Repetition Frequency (Hz)
pulseWidthUp = 100e-6;      % Up-chirp duration (sec)
startFreqUp = 3e4;          % Up-chirp start freq (Hz)
stopFreqUp = 6e4;           % Up-chirp end freq (Hz)
pulseWidthDown = 200e-6;    % Down-chirp duration (sec)
startFreqDown = 6e4;        % Down-chirp start freq (Hz)
stopFreqDown = 3e4;         % Down-chirp end freq (Hz)

%% Scenario Parameters
numPulses = 10;           % Number of consecutive pulse segments to sim
maxEchoesPerSegment = 1;    % Maximum number of echoes per segment
snrDb = -10;                % Received SNR of echoes
velocityMph = 5000;         % Radial velocity (Mph), positive inward

%% Generate simulation objects
cfgObj = Config(sampleRate, prf, pulseWidthUp, startFreqUp, ...
    stopFreqUp, pulseWidthDown, startFreqDown, stopFreqDown);
txRxObj = Transceiver(cfgObj, numPulses=numPulses, matchFiltType="fft");
sceneObj = Scenario(cfgObj, txRxObj, maxEchoesPerSegment, ...
    numPulses=numPulses, snrDb=snrDb, velocityMph=velocityMph);

%% Simulate
txSamps = txRxObj.generate(); % generates transmitted samples
rxSamps = sceneObj.apply();    % creates received samples with echos
mfOut = txRxObj.matchFilter(rxSamps); % applies matched filter and segments

%% Plots
Visualizer.chirps(txRxObj);
Visualizer.segments(cfgObj, sceneObj, txRxObj, rxSamps, mfOut);