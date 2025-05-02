classdef Scenario < handle
    % This class generates channel scenarios based on user defined input

    properties
        cfgObj Config
        txRxObj Transceiver
        randomEchos (1,1) logical = true
        userEchoDelays (:,1) {mustBeNonnegative} = []
        echoDelays (:,:) {mustBeFloat}
        echoTaps (:,:) {mustBeFloat}
        numPulses (1,1) {mustBeFloat}
        snrDb (:,1) {mustBeFloat}
        minSnrDb (1,1) {mustBeFloat} = 0
        maxSnrDb (1,1) {mustBeFloat} = 15
        maxEchoesPerSegment (1,1) {mustBeFloat}
        dopplerFreq (1,1) {mustBeFloat}
    end
    methods
        function this = Scenario(cfgObj, txRxObj, maxEchoesPerSegment, opts)
            arguments(Input)
                cfgObj Config
                txRxObj Transceiver
                maxEchoesPerSegment (1,1) {mustBeFloat} = 1
                opts.numPulses (1,1) {mustBeFloat} = 10
                opts.snrDb (:,1) {Scenario.mustBeCorrectSize(opts.snrDb, ...
                    maxEchoesPerSegment)} = []
                opts.velocityMph (1,1) {mustBeFloat} = 0
                opts.echoDelays (:,1) {mustBeNonnegative} = []
            end
            this.cfgObj = cfgObj;
            this.txRxObj = txRxObj;
            this.maxEchoesPerSegment = maxEchoesPerSegment;
            this.numPulses = opts.numPulses;
            if isempty(opts.echoDelays)
                this.randomEchos = true;
            else
                this.randomEchos = false;
                this.userEchoDelays = opts.echoDelays;
            end
            if isempty(opts.snrDb)
                this.snrDb = (this.maxSnrDb+this.minSnrDb) * ...
                    rand(maxEchoesPerSegment,1) - this.minSnrDb;
            else
                this.snrDb = opts.snrDb;
            end
            velocity = opts.velocityMph * 0.44704;
            this.dopplerFreq = 2 * velocity * ...
                this.cfgObj.centerFreq / this.cfgObj.lightspeed;
        end

        function r = apply(this)
            % Inject echoes for each PRI with FIR model

            if this.randomEchos
                this.echoDelays = randi([1, this.cfgObj.maxPulseDelay], ...
                    this.numPulses, this.maxEchoesPerSegment);
            else
                this.echoDelays = this.userEchoDelays.' .* ...
                    ones(this.numPulses, this.maxEchoesPerSegment);
            end

            this.echoTaps = sqrt(10.^(this.snrDb'/10)) .* ...
                exp(1j*2*pi*rand(size(this.echoDelays,1),1));

            r = zeros(this.numPulses*this.cfgObj.sampsPerSegment,1);
            for k = 1:this.numPulses
                startIdx = (k-1)*this.cfgObj.sampsPerSegment + 1;
                d = this.echoDelays(k,:);
                if mod(k,2)==1
                    pulse = this.txRxObj.upChirp;
                else
                    pulse = this.txRxObj.downChirp;
                end
                channFilt = zeros(this.cfgObj.sampsPerSegment,1);
                channFilt(d) = this.echoTaps(k,:);
                seg = conv(channFilt, pulse, 'same');
                r(startIdx:startIdx+this.cfgObj.sampsPerSegment-1) = seg;
            end

            r = r .* exp(1j*2*pi* ...
                (this.dopplerFreq/this.cfgObj.sampleRate) ...
                *(0:length(r)-1).');
            r = r + sqrt(1/2)*(randn(size(r))+1j*randn(size(r)));
        end
    end

    methods(Static, Access = private)
        function mustBeCorrectSize(in, ref)
            if ~( isscalar(in) || isequal(numel(in), ref) || isempty(in))
                error( ...
                    'MyClass:Scenario', ...
                    'snrDb must be a scalar or the same size as maxEchoes.' ...
                    );
            end
        end
    end
end