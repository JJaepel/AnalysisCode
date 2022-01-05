function driftingGrating_Fct(analysisParams)
close all

%% 0.) define folders and structures
analysisParams = createFolderEpi(analysisParams);
%clear out all files in folder
if analysisParams.clean
    disp('Deleting old analysis files')
    delete([analysisParams.saveDirectory '\*'])
end

%% 1.) Load tiffs
t0=tic;
data.rawF = readingImagingData(analysisParams.EpiDirectory, analysisParams.downsample);
toc(t0)

%% 2.) load metadata
metadata.StimParams=LoadStimParams(analysisParams.Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(analysisParams.Sp2dDirectory);
metadata.StimParams.path=fullfile(analysisParams.Sp2dDirectory);

%create StimCodes
metadata = createStimCodesEpi(metadata);

%% 3.) get ROI of window and of subregions
analysis = struct;
analysis.rawFMeanImg = mean(data.rawF,3);
analysis.baseImg = mean(data.rawF(:,:,1:50),3);
analysis.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);
analysis.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
analysis = makeMasks(data, analysis, analysisParams,1);

%% 4.) Chop traces
disp('Chopping Traces')
[analysis, metadata] = ChopStimulusTraceEpi(analysis,metadata,data,analysisParams.field, analysisParams.bandpass);
clear data

%% 5.) Make timecourse
disp('Plotting time course')
makeTimecourseEpi(metadata, analysis, analysisParams)

%% 6.) Map analysis
analysis = showEpiRespAvg(analysis, metadata, analysisParams.field, analysisParams.saveDirectory);

%% 7.) Add contours for single trial analysis
analysis = makeContourMaps(metadata, analysis, analysisParams);

%% 8.) Save data
disp('Saving analyzed data')
save(fullfile(analysisParams.saveDirectory, 'AnaData.mat'),'-v7.3', 'metadata', 'analysisParams', 'analysis');


