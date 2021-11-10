function spontaneous_Fct(analysisParams)
close all

%% 0.) define folders and structures
analysisParams = createFolderEpi(analysisParams);
%clear out all files in folder
if analysisParams.clean
    disp('Deleting old analysis files')
    delete([analysisParams.saveDirectory '\*'])
end

%% 1.) load metadata
try
    metadata.Imaging=LoadFrameTimes(analysisParams.Sp2dDirectory);
catch
    metadata.Imaging.rate = 14.99;
end


%% 2.) Load tiffs
disp('Loading raw data')
data.rawF = readingImagingData(analysisParams.EpiDirectory);

% if applicaple, downsample
if analysisParams.downsample > 1
    data = downSampleData(data, analysisParams.downsample);
end

% if applicable, register data to remove movement
if analysisParams.register
    numberOfImages = size(data.rawF,3);
    try
        disp('Loading shift data')
        load(fullfile(analysisParams.EpiDirectory, 'shifts.mat'),'rowShift', 'colShift')
    catch
        referenceImg = real(bandpassFermiFilter_Revised(data.rawF(:,:,1),-1,600,1000/172));
        rowShift = zeros(numberOfImages,1);
        colShift = zeros(numberOfImages,1);
        tic
        disp('Registering images')
        parfor(ii = 2:numberOfImages) %takes roughly 15 min
            [~, rowShift(ii), colShift(ii)] = registerImages(referenceImg,data.rawF(:,:,ii),[-1 600 1000/172],true)
        end
        toc
        save(fullfile(analysisParams.EpiDirectory, 'shifts.mat'),'rowShift', 'colShift');
    end
    disp('Shifting images')
    for ii = 2:numberOfImages
        [data.rawF(:,:,ii)] = shiftImages(data.rawF(:,:,ii),rowShift(ii), colShift(ii), true);
    end
    view_tiff_stack(data.rawF, metadata)
end

%% 3.) get ROI of window and of subregions
analysis = struct;
analysis.rawFMeanImg = mean(data.rawF,3);
analysis.baseImg = mean(data.rawF(:,:,1:50),3);
analysis.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);
analysis.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 

disp('Making masks')
analysis = makeMasks(data, analysis, analysisParams);

%% 4.) Apply dff
disp('Normalizing Stack')
data = NormalizeStack(metadata,data); %takes roughly 10 min
data.dff(~analysis.maskBV(:)) = NaN;

disp('Z-scoring Stack')
mu = nanmean(data.dff(:));
sd = nanstd(data.dff(:));
data.zscore = (data.dff - mu)/sd;
data.zscorePCA = data.zscore;
data.zscorePCA(isnan(data.zscore)) = 0;

[~, mixedfilters, ~] = SVDsimple(data.zscorePCA,1);
data.PCAsZScore = mixedfilters;

figure
imagesc(data.PCAsZScore(:,:,1))
saveas(gcf, [analysisParams.saveDirectory, 'Mask1stPCA_ZScore.png'])
close gcf

%% 5.) Get active frames
disp('Detecting active frames')
[data.activeFrames, analysis.maskBVActive, analysis.eventOnsets] = getActiveFramesMasked(data.dff,analysis.maskBV);

%% 6.) If applicable, apply bandpassfilter
if analysisParams.bandpass
    disp('Applying bandpass filter')
    analysis.ROIActive =true( [size(data.activeFrames,1),size(data.activeFrames,2)]);
    data.activeFrames = LowHighNormalize(double(data.activeFrames), analysis.ROIActive);
end

%% 7.) get ative Frames, compute correlations of the imaging stack and show it
% Computes correlations of the imaging stack (spontaneous, response, signal, or noise). 
            % The image dimensions can be "N" dimensions, but must be organized such that they are:
            %    *Spontaneous: Should be (x,y,t). The active frames are automatically extracted. 
            %    *Response:    Collapsed from a multidimensional array (i.e. [x,y,nCond,nTrials,t])
            %                  to (x,y,n), where n is all images
            %    *Signal:      Takes in a (x,y,nCond,nTrials,t), averages along the fourth dimension,
            %    *Noise        Computes correlations along nCond.

analysis.corrTable = computeCorrelationTable(data.activeFrames,analysis.ROIActive);
%showCorrelationStructure(analysis.corrTable,analysis, analysisParams.saveDirectory)

%% 8.) Calculate fracture lines
%temporaly clear workspace
corrTable = analysis.corrTable;
activeFrames = data.activeFrames(:,:,1);
ROIActive = analysis.ROIActive;
save(fullfile(analysisParams.saveDirectory, 'AnaData.mat'),'-v7.3', 'metadata', 'analysisParams', 'analysis');

clear data
clear analysis
clear mixedfilters

disp('Calculating fracture maps')
analysisParams.sizeMap = 3; 
fractureMap = makeFractureLines(ROIActive, activeFrames, corrTable, analysisParams.sizeMap, analysisParams);

%% 9.) Save data
disp('Saving analyzed data')
save(fullfile(analysisParams.saveDirectory, 'AnaData.mat'),'-v7.3', 'metadata', 'analysisParams', 'analysis');

