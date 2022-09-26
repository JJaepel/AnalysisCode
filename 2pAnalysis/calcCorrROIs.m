function [corrROIs] = calcCorrROIs(analysisParams, analysis)

%make empty array 
corrROIs= zeros(length(analysis.(analysisParams.field).roi), length(analysis.(analysisParams.field).roi));

%load traces for each ROI and calculate the correlation between the traces
for A = 1:length(analysis.(analysisParams.field).roi)
    traceA = analysis.(analysisParams.field).roi(A).normResponse;
    for B = 1:length(analysis.(analysisParams.field).roi)
        traceB = analysis.(analysisParams.field).roi(B).normResponse;
        corrROIs(A,B) = corr(traceA', traceB');
    end
end
