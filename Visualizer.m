classdef Visualizer < handle
    % This class holds plotting functions that are useful for LFM radar
    
    methods(Static)
        function chirps(txRxObj)
            t = tiledlayout(2, 1, 'TileSpacing', 'compact', ...
                'Padding', 'compact');

            nexttile
            plot(real(txRxObj.upChirp)); hold all
            plot(imag(txRxObj.upChirp))
            title('Up Chirp');
            legend('Real', 'Imag')

            nexttile
            plot(real(txRxObj.downChirp)); hold all
            plot(imag(txRxObj.downChirp))
            title('Down Chirp');
            legend('Real', 'Imag')

            xlabel(t, 'Time (s)')
            ylabel(t, 'Amplitude')
        end

        function segments(cfgObj, sceneObj, txRxObj, samples, mfSegments)
            sampsPerSegment = cfgObj.sampsPerSegment;
            segLines = (0:sampsPerSegment:length(samples)).';
            delays = sceneObj.echoDelays + segLines(1:end-1);
            thresholds = txRxObj.getThresholds(1e-6);

            figure;
            subplot(2,1,1)
            plot(real(samples)); hold on            
            for k = 1:length(segLines)
                xline(segLines(k), '--', sprintf('Segment %i', k), ...
                    'Color',[.4 .4 .4], 'LabelVerticalAlignment','top');
            end
            title('Matched Filter Input (Real)');
            ylabel('Amplitude')
            xlabel('Sample Delay')
            legend('MF Out', 'location', 'southeast')

            subplot(2,1,2)
            h1 = semilogy(abs(mfSegments(:)).^2); hold on
            h2 = plot(delays(:), abs(mfSegments(delays(:))).^2, 'o');
            axis([-inf inf 10^3 10^10])
            h4 = [];
            for k = 1:length(segLines)-1                
                x = segLines(k):segLines(k+1);
                if mod(k,2) == 1
                    h3 = plot(x, thresholds(1)*ones(size(x)), ...
                        'k-.', 'LineWidth', 1 );
                    xline(segLines(k), '--', sprintf('Segment %i\nUp    ', k), ...
                    'Color',[.4 .4 .4], 'LabelVerticalAlignment','top');
                else
                    h4 = plot(x, thresholds(2)*ones(size(x)), ...
                        'k:', 'LineWidth', 1 );
                    xline(segLines(k), '--', sprintf('Segment %i\nDown  ', k), ...
                    'Color',[.4 .4 .4], 'LabelVerticalAlignment','top');
                end
            end
            title('Matched Filter Output');
            xlabel('Sample Delay')
            ylabel('|MF|^2')
            if isempty(h4)
                legend([h1 h2 h3], 'MF Out', 'True Delays', ...
                    'Up Thresh', 'location', 'southeast')
            else
                legend([h1 h2 h3 h4], 'MF Out', 'True Delays', ...
                    'Up Thresh', 'Down Thresh', 'location', 'southeast')
            end
        end
    end
end
