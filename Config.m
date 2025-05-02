classdef Config < handle
    % This class holds the user defined configuration parameters that
    % define the LFM radar properties. It allows us to validate input early
    % to catch poorly conceived parameter selections. It also allows
    % derived quantities to be computed and stored.

    properties
        sampleRate (1,1) {mustBeFloat}     % Sampling rate (Hz)
        prf (1,1) {mustBeFloat}            % Pulse repetition frequency (Hz)
        pri (1,1) {mustBeFloat}            % Pulse repetition interval (sec)
        pulseWidthUp (1,1) {mustBeFloat}   % Up-chirp duration (s)
        startFreqUp (1,1) {mustBeFloat}    % Up-chirp start frequency (Hz)
        stopFreqUp (1,1) {mustBeFloat}     % Up-chirp end frequency (Hz)
        pulseWidthDown (1,1) {mustBeFloat} % Down-chirp duration (s)
        startFreqDown (1,1) {mustBeFloat}  % Down-chirp start frequency (Hz)
        stopFreqDown (1,1) {mustBeFloat}   % Down-chirp end frequency (Hz)
        sampsPerSegment (1,1) {mustBeFloat}% Samples per segment
        maxPulseDelay (1,1) {mustBeFloat}  % Max delay of any pulse
        centerFreq (1,1) {mustBeFloat}   % Center frequency of chirps
        lightspeed (1,1) {mustBeFloat} = 299792458
    end
    methods
        function this = Config(sampleRate, prf, pulseWidthUp, ...
                startFreqUp, stopFreqUp, pulseWidthDown, ...
                startFreqDown, stopFreqDown)
            arguments(Input)
                sampleRate (1,1) {mustBeFloat}
                prf (1,1) {mustBeFloat}
                pulseWidthUp (1,1) {mustBeFloat}
                startFreqUp (1,1) {mustBeFloat}
                stopFreqUp (1,1) {mustBeFloat}
                pulseWidthDown (1,1) {mustBeFloat}
                startFreqDown (1,1) {mustBeFloat}
                stopFreqDown (1,1) {mustBeFloat}
            end
            this.sampleRate = sampleRate;
            this.prf = prf;
            this.pulseWidthUp = pulseWidthUp;
            this.startFreqUp = startFreqUp;
            this.stopFreqUp = stopFreqUp;
            this.pulseWidthDown = pulseWidthDown;
            this.startFreqDown = startFreqDown;
            this.stopFreqDown = stopFreqDown;            

            % Derived quantities
            this.pri = 1/this.prf;
            this.sampsPerSegment = round(sampleRate / prf);
            this.maxPulseDelay = (this.pri - ...
                max(pulseWidthUp, pulseWidthDown)) * sampleRate;
            this.centerFreq = (this.startFreqDown + this.stopFreqDown)/2;
        end
    end
end