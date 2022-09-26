function analysis = classifyRespROIs(analysisParams, analysis, metadata, ROI)

pretrialTime= analysis.(analysisParams.field).preTrialTime;
preTrialIndex= (1:floor(pretrialTime * metadata.TwoPhoton.rate));
baselines=analysis.(analysisParams.field).roi(ROI).stimResponseTrace(:,:,preTrialIndex);
analysis.(analysisParams.field).roi(ROI).baselineSD = std(baselines,[],3);
analysis.(analysisParams.field).roi(ROI).baselineMean = mean(baselines,3);
analysis.(analysisParams.field).roi(ROI).baselineMean(analysis.(analysisParams.field).roi(ROI).baselineMean < 0) = mean(mean(analysis.(analysisParams.field).roi(ROI).baselineMean,2),1);
analysisPeriod=(analysis.(analysisParams.field).windowStart:analysis.(analysisParams.field).windowStop);
stimResp = analysis.(analysisParams.field).roi(ROI).stimResponseTrace(:,:,analysisPeriod);
analysis.(analysisParams.field).roi(ROI).peaks = max(stimResp,[],3);
analysis.(analysisParams.field).roi(ROI).zscore = ([analysis.(analysisParams.field).roi(ROI).peaks]-[analysis.(analysisParams.field).roi(ROI).baselineMean])./[analysis.(analysisParams.field).roi(ROI).baselineSD];
analysis.(analysisParams.field).roi(ROI).crosser = sum(analysis.(analysisParams.field).roi(ROI).zscore > analysisParams.zThresh,2);
analysis.(analysisParams.field).roi(ROI).respStim = analysis.(analysisParams.field).roi(ROI).crosser >= ((metadata.StimParams.numTrials)*analysisParams.fraction);
if sum(analysis.(analysisParams.field).roi(ROI).respStim) > 0
    analysis.(analysisParams.field).roi(ROI).isResponseSignificant = 1;
else 
    analysis.(analysisParams.field).roi(ROI).isResponseSignificant = 0;
end

if analysisParams.stimType == 2 || analysisParams.stimType == 3 || analysisParams.stimType == 5
    for con = 1:metadata.StimParams.numCon
        tempCrosser = analysis.(analysisParams.field).roi(ROI).crosser((con-1)*metadata.StimParams.numOrientations+1:con*metadata.StimParams.numOrientations);
        tempRespStim = tempCrosser > ((metadata.StimParams.numTrials)*analysisParams.fraction);
        if sum(tempRespStim) > 0
            isResp_allCon(con) = 1;
        else
            isResp_allCon(con) = 0;
        end
    end
    analysis.(analysisParams.field).roi(ROI).isResp_allCon = isResp_allCon;

end