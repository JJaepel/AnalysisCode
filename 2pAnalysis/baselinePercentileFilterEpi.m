function highpassFilteredTrace = baselinePercentileFilterEpi(inputTrace, fps, filteredCutoff, desiredPercentileRank)            

    if nargin < 3
        filteredCutoff = 60;
    end
    if nargin < 4
        desiredPercentileRank = 20;
    end

    paddingLength = ceil(length(inputTrace)/1);
    paddedTrace   = [inputTrace(paddingLength:-1:1); inputTrace; inputTrace(paddingLength:-1:1)];
    filteredTrace = percentileFilt1(paddedTrace, desiredPercentileRank, round(filteredCutoff*fps));
    filteredTrace = filteredTrace(paddingLength+(1:length(inputTrace)));
    
    % The low-pass filter the filtered trace to smooth it out
    butterWorthOrder = 1;
    Wn = (1/filteredCutoff) / (fps/2);
    [b,a] = butter(butterWorthOrder, Wn, 'low');
    highpassFilteredTrace = filtfilt(b,a,[filteredTrace(paddingLength:-1:1); filteredTrace; filteredTrace(1:paddingLength)]);
    highpassFilteredTrace = highpassFilteredTrace(paddingLength+[1:length(inputTrace)]);
end