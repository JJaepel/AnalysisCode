function HLSmap = stimulusMap(stimulusAveragedMaps, clipValue)
    % Calculate max stimulus, amplitude and tuning width (circular variance of sorted map)
    NumStims = size(stimulusAveragedMaps, 3);
   [ResponseAmplitude,PrefStim] = max(stimulusAveragedMaps, [], 3 );

    % Create HLS map
    if NumStims > 3
        figure;
        CM = colormap('jet');
        close;
        IX = round(linspace(1,64,NumStims));
        for i = 1:length(IX)
            C{i} = CM(IX(i),:);
        end
    elseif NumStims == 3
        C{1} = [1 0 0];
        C{2} = [0 1 0];
        C{3} = [0 0 1];
    elseif NumStims == 2
        C{1} = [0 0.8 1];
        C{2} = [1 0 0.4];
    end
    
    excludingNaNs = ~isnan(ResponseAmplitude(:));
    highClipVal = mean(ResponseAmplitude(excludingNaNs(:)))+clipValue(2)*std(ResponseAmplitude(excludingNaNs(:)));
    lowClipVal  = mean(ResponseAmplitude(excludingNaNs(:)))-clipValue(1)*std(ResponseAmplitude(excludingNaNs(:)));
    ResponseAmplitude(ResponseAmplitude > highClipVal) = highClipVal;
    ResponseAmplitude(ResponseAmplitude < lowClipVal ) = lowClipVal;
    ResponseAmplitude_Scaled = rescale(ResponseAmplitude);
    
    HLSmap = zeros( size(stimulusAveragedMaps,1), size(stimulusAveragedMaps,2), 3 );            
    for y = 1:size(stimulusAveragedMaps,1)
        for x = 1:size(stimulusAveragedMaps,2)
            for c = 1:3
                if ~isnan(PrefStim(y,x)) && ResponseAmplitude_Scaled(y,x) > 0
                    HLSmap(y,x,c) = (C{PrefStim(y,x)}(c) .* ResponseAmplitude_Scaled(y,x));
                end
            end
        end
    end

    % Add color index
    [yRes,xRes,~] = size(HLSmap);
    xRes = xRes - 100;
    for s = 1:NumStims
        xRange = round(50+(((s*(xRes/NumStims)) - (xRes/(NumStims*2))):(s*(xRes/NumStims))));
        for c = 1:3
            HLSmap(yRes:yRes+10,xRange,c) = C{s}(c);
        end
    end
end