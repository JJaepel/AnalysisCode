function spontaneous_Fct(analysisParams)
close all

%% 0.) define folders and structures
analysisParams = createFolderEpi(analysisParams);
analysisParams.clean = 0;
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
rawF = readingImagingData(analysisParams.EpiDirectory);

% if applicaple, downsample
if analysisParams.downsample > 1
    rawF = downSampleData(rawF, analysisParams.downsample);
end

% if applicable, register data to remove movement
if analysisParams.register
    numberOfImages = size(rawF,3);
    try
        disp('Loading shift data')
        load(fullfile(analysisParams.EpiDirectory, 'shifts.mat'),'rowShift', 'colShift')
    catch
        referenceImg = real(bandpassFermiFilter_Revised(rawF(:,:,1),-1,600,1000/172));
        rowShift = zeros(numberOfImages,1);
        colShift = zeros(numberOfImages,1);
        tic
        disp('Registering images')
        parfor(ii = 2:numberOfImages) %takes roughly 15 min
            [~, rowShift(ii), colShift(ii)] = registerImages(referenceImg,rawF(:,:,ii),[-1 600 1000/172],true)
        end
        toc
        save(fullfile(analysisParams.EpiDirectory, 'shifts.mat'),'rowShift', 'colShift');
    end
    disp('Shifting images')
    for ii = 2:numberOfImages
        [rawF(:,:,ii)] = shiftImages(rawF(:,:,ii),rowShift(ii), colShift(ii), true);
    end
    view_tiff_stack(rawF, metadata)
end

%% 3.) get ROI of window and of subregions
analysis = struct;
analysis.rawFMeanImg = mean(rawF,3);
analysis.baseImg = mean(rawF(:,:,1:50),3);
analysis.gaussMeanImg = imgaussfilt(mean(rawF, 3), 4);
analysis.ROI =true( [size(rawF,1),size(rawF,2)]); 

disp('Making masks')
analysis = makeMasks(rawF, analysis, analysisParams.saveDirectory,0);

%% 4.) Apply dff
disp('Normalizing Stack')
data.dff = NormalizeStack(metadata,rawF); %takes roughly 10 min
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
[data.activeFrames, analysis.maskBVActive, analysis.eventOnsets] = getActiveFramesMasked(data.zscore,analysis.maskBV);

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

analysis.corrTable = computeCorrelationTable(data.activeFrames,analysis.maskBVActive);
showCorrelationStructure(analysis.corrTable,analysis.maskBVActive, analysis, analysisParams.saveDirectory)

%% 8.) Calculate fracture lines and saving data
%temporaly clear workspace
corrTable = analysis.corrTable;
activeFrames = data.activeFrames;
ROIActive = analysis.ROIActive;
save(fullfile(analysisParams.saveDirectory, 'AnaData.mat'),'-v7.3', 'metadata', 'analysisParams', 'analysis');

clear data
clear analysis
clear mixedfilters

disp('Calculating fracture maps')
sizeOptions = [1,3,5,7];
for sz = 1:length(sizeOptions)
    fractureMap{sz} = makeFractureLines(ROIActive, activeFrames, corrTable, sizeOptions(sz), analysisParams);
end

%% 9.) Save data
disp('Saving analyzed data')
save(fullfile(analysisParams.saveDirectory, 'Fractures.mat'),'-v7.3', 'fractureMap');

