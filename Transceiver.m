classdef Transceiver < handle
    properties
        cfgObj Config
        sceneObj Scenario
        upChirp (:,1) {mustBeFloat}
        upChirpMf (:,1) {mustBeFloat}
        numSampsUpChirp (1,1) {mustBePositive} = 1
        downChirp (:,1) {mustBeFloat}
        downChirpMf (:,1) {mustBeFloat}
        numSampsDownChirp (1,1) {mustBePositive} = 1
        numPulses (1,1) {mustBePositive} = 1

        % Strategy pattern: handle to selected match filter function
        matchFilterFcn function_handle

        % Precomputed FFT resources
        nfft (1,1) {mustBePositive} = 1
        upChirpMfFreq (:,1) {mustBeFloat}
        downChirpMfFreq (:,1) {mustBeFloat}
        startIdxUp (1,1) {mustBePositive} = 1
        startIdxDown (1,1) {mustBePositive} = 1

        Hmat (:,:) {mustBeFloat}
        oddIdxs (:,1) {mustBePositive} = 1
        evenIdxs (:,1) {mustBePositive} = 1

        % Preallocated matched filter memory
        mfOut (:,:) {mustBeFloat}
        segments (:,:) {mustBeFloat}
    end

    methods
        function this = Transceiver(cfgObj, opts)
            arguments(Input)
                cfgObj Config
                opts.matchFiltType string ...
                    {Transceiver.mustBeMfType(opts.matchFiltType)} = "fft"
                opts.numPulses (1,1) {mustBePositive} = 10
            end
            this.cfgObj = cfgObj;
            this.numPulses = opts.numPulses;

            % Generate up-chirp and its matched filter
            this.upChirp = Transceiver.makeChirp(cfgObj.sampleRate, ...
                cfgObj.pulseWidthUp, cfgObj.startFreqUp, ...
                cfgObj.stopFreqUp);
            this.numSampsUpChirp = length(this.upChirp);
            this.upChirpMf = Transceiver.makeMatchedFilter(this.upChirp);

            % Generate down-chirp and its matched filter
            this.downChirp = Transceiver.makeChirp(cfgObj.sampleRate, ...
                cfgObj.pulseWidthDown, cfgObj.startFreqDown, ...
                cfgObj.stopFreqDown);
            this.numSampsDownChirp = length(this.downChirp);
            this.downChirpMf = Transceiver.makeMatchedFilter(this.downChirp);

            % Preallocate for efficiency
            this.mfOut = zeros(cfgObj.sampsPerSegment, this.numPulses);
            this.segments = zeros(this.cfgObj.sampsPerSegment, ...
                    this.numPulses);

            % Precompute as much matched filter values as possible
            this.oddIdxs  = 1:2:this.numPulses;
            this.evenIdxs = 2:2:this.numPulses;
            if strcmpi(opts.matchFiltType, "fft")
                maxLen = max(this.numSampsUpChirp, this.numSampsDownChirp);
                this.nfft = 2^nextpow2(cfgObj.sampsPerSegment + maxLen - 1);
                this.upChirpMfFreq = fft(this.upChirpMf, this.nfft);
                this.downChirpMfFreq = fft(this.downChirpMf, this.nfft);
                this.startIdxUp = floor(this.numSampsUpChirp/2) + 1;
                this.startIdxDown = floor(this.numSampsDownChirp/2) + 1;
                
                this.Hmat = zeros(this.nfft, this.numPulses);
                this.Hmat(:,this.oddIdxs) = ...
                    repmat(this.upChirpMfFreq, 1, numel(this.oddIdxs));
                this.Hmat(:,this.evenIdxs) = ...
                    repmat(this.downChirpMfFreq, 1, numel(this.evenIdxs));

                this.matchFilterFcn = @(in) this.matchFilterFftImpl(in);
            else
                this.matchFilterFcn = @(in) this.matchFilterConvImpl(in);
            end
        end

        function tx = generate(this)
            % Create alternating up/down chirp train
            
            tx = zeros(this.numPulses * this.cfgObj.sampsPerSegment, 1);
            for k = 1:this.numPulses
                startIdx = (k-1)*this.cfgObj.sampsPerSegment + 1;
                if mod(k,2)==1
                    wave = this.upChirp;
                else
                    wave = this.downChirp;
                end
                tx(startIdx:startIdx+length(wave)-1) = wave;
            end

            tx = tx/sqrt(mean(abs(tx).^2));
        end

        function out = matchFilter(this, samples)
            % Apply matched filtering segment-by-segment. Using the
            % strategy pattern allows us to present a unified method call
            % to the user based on their selection during object creation
            arguments(Input)
                this
                samples (:,1) {mustBeFloat}
            end

            % Delegate to the chosen implementation
            out = this.matchFilterFcn(samples);
        end

        function thresholds = getThresholds(this, pFa, noisePow)
            % Thresholds that will achieve target prob false alarm
            % The first element is the upChirp threshold, second is down
            arguments(Input)
                this
                pFa (1,1) {mustBeFloat}
                noisePow (1,1) {mustBeFloat} = 1
            end
            thresholds(1) = noisePow * -log(pFa) * ...
                sum(abs(this.upChirpMf).^2);
            thresholds(2) = noisePow * -log(pFa) * ...
                sum(abs(this.downChirpMf).^2);
        end
    end

    methods(Access=private)

        function out = matchFilterFftImpl(this, samples)
            % Apply matched filtering segment-by-segment

            this.segments(:) = samples; % use preallocated memory            
            X = fft(this.segments, this.nfft, 1);
            Y = X .* this.Hmat;
            yFull = ifft(Y, this.nfft, 1);

            this.mfOut(:,this.oddIdxs) = ...
                yFull(this.startIdxUp:this.startIdxUp + ...
                this.cfgObj.sampsPerSegment-1, this.oddIdxs);
            this.mfOut(:,this.evenIdxs) = ...
                yFull(this.startIdxDown:this.startIdxDown + ...
                this.cfgObj.sampsPerSegment-1, this.evenIdxs);
            out = this.mfOut;
        end

        function out = matchFilterConvImpl(this, samples)
            % Apply matched filtering segment-by-segment

            this.segments(:) = samples;

            this.mfOut(:,this.oddIdxs) = ...
                conv2(this.segments(:,this.oddIdxs), ...
                this.upChirpMf, 'same');

            this.mfOut(:,this.evenIdxs) = ...
                conv2(this.segments(:,this.evenIdxs), ...
                this.downChirpMf, 'same');

            out = this.mfOut;
        end
    end

    methods (Static, Access = private)
        function out = makeChirp(sampleRate, pulseWidth, ...
                startFreq, stopFreq)
            % Create a chirp with user defined properties
            time = (0 : 1/sampleRate : pulseWidth-1/sampleRate).';
            freqSlope = (stopFreq - startFreq) / pulseWidth;
            out = exp(1j*2*pi*(startFreq*time + ...
                0.5*freqSlope*time.^2));            
        end

        function out = makeMatchedFilter(in)
            % Generate a matched filter from the input
            out = conj(flipud(in));
        end

        function mustBeMfType(x)
            % Custom validator function for user input check
            if ~any(strcmpi(x, ["fft","conv"]))
                error( ...
                  'Transceiver:InvalidMatchFiltType', ...
                  "matchFiltType must be 'fft' or 'conv' (case-insensitive)" ...
                );
            end
        end
    end
end