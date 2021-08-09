close all; 
clear all
imageDirectory = 'Z:\Juliane\Data\Epi\'; 
metaDirectory = 'Z:\Juliane\Data\Spike2Data\';
name = 'F2500_2021-03-17'; 
exptID = 2;
sp2ID = 2;
spatialDownsamplingFactor = 1;
analysisDir = 'Z:\Juliane\Data\ImageAnalysis\';
timecourse = 1;
responseSign = 1;
contours = 0;

%define location of imaging data and meta data
imagingDataDirectory = sprintf('%s%s\\tseries_%d\\',imageDirectory,name,exptID);
Sp2dDirectory    = sprintf('%s%s\\t%05d\\',metaDirectory,name,sp2ID);
saveDataDirectory   = sprintf('%s%s\\t%05d\\',analysisDir,name,exptID);
saveDataDir = sprintf('%s%s\\',analysisDir,name);

%make savedir
if ~exist(saveDataDirectory, 'dir')
    mkdir(saveDataDirectory);
end

%% load imaging data
tifStack = LoadImagingData(imagingDataDirectory, spatialDownsamplingFactor);
expParam.ROI =true( [size(tifStack,1),size(tifStack,2)]); 
expParam.rawFMeanImg = mean(tifStack,3);
expParam.baseImg = mean(tifStack(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(tifStack, 3), 4);

%% load metadata
metadata.StimParams=Load_stimparams(Sp2dDirectory);
metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
metadata.StimParams.path=fullfile(Sp2dDirectory);
metadata.StimParams.series=expt_id;
switch metadata.StimParams.type
    case 'continuousEdge_withOrientationsAndTriggers'
        if isa(metadata.StimParams.startPoint, 'double')
            if metadata.StimParams.orientations == 90 || metadata.StimParams.orientations == 270
                metadata.StimParams.visualSpace = 'azimuth';
            elseif metadata.StimParams.orientations == 0 || metadata.StimParams.orientations == 180
                metadata.StimParams.visualSpace = 'elevation';
            end
        else
            if isempty(metadata.StimParams.right_screen)
                metadata.StimParams.startPoint = str2double(regexp({metadata.StimParams.startPoint},'-?[0-9]*', 'match', 'once'));
                metadata.StimParams.endPoint = str2double(regexp({metadata.StimParams.endPoint},'-?[0-9]*', 'match', 'once'));
            else
                if metadata.StimParams.orientations == 90
                    metadata.StimParams.startPoint = -metadata.StimParams.left_screen;
                    metadata.StimParams.endPoint = metadata.StimParams.right_screen;
                    metadata.StimParams.visualSpace = 'azimuth';
                elseif metadata.StimParams.orientations == 270
                    metadata.StimParams.startPoint = metadata.StimParams.right_screen;
                    metadata.StimParams.endPoint = -metadata.StimParams.left_screen;
                    metadata.StimParams.visualSpace = 'azimuth';
                elseif metadata.StimParams.orientations == 0
                    metadata.StimParams.startPoint = -metadata.StimParams.lower_screen;
                    metadata.StimParams.endPoint = metadata.StimParams.upper_screen; 
                    metadata.StimParams.visualSpace = 'elevation';
                elseif metadata.StimParams.orientations == 180
                    metadata.StimParams.startPoint = metadata.StimParams.upper_screen;
                    metadata.StimParams.endPoint =-metadata.StimParams.lower_screen; 
                    metadata.StimParams.visualSpace = 'elevation';
                end
            end
        end
    case 'continuousBar_withTriggers'
        if metadata.StimParams.stimID == 1
            metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
            metadata.StimParams.endPoint = - metadata.StimParams.endPoint; 
            metadata.StimParams.visualSpace = 'elevation';
        elseif metadata.StimParams.stimID == 3
            metadata.StimParams.visualSpace = 'azimuth';
        elseif    metadata.StimParams.stimID == 5
            metadata.StimParams.visualSpace = 'elevation';
        elseif   metadata.StimParams.stimID == 7
            metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
            metadata.StimParams.endPoint = - metadata.StimParams.endPoint; 
            metadata.StimParams.visualSpace = 'azimuth';
        end        
    case 'continuousBar_withTriggers_down'
        metadata.StimParams.orientations = 180;
        metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
        metadata.StimParams.endPoint = - metadata.StimParams.endPoint;
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousBar_withTriggers_up'
        metadata.StimParams.orientations = 0;
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousBar_withTriggers_left'
        metadata.StimParams.orientations = 90;
        metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
        metadata.StimParams.endPoint = - metadata.StimParams.endPoint;
        metadata.StimParams.visualSpace = 'azimuth';
    case 'continuousBar_withTriggers_right'
        metadata.StimParams.orientations = 270;
        metadata.StimParams.visualSpace = 'azimtuh';
    case 'continuousEdge_azimuth'
        if metadata.StimParams.orientations == 90
           metadata.StimParams.startPoint = -metadata.StimParams.left_screen;
           metadata.StimParams.endPoint = metadata.StimParams.right_screen;
        elseif metadata.StimParams.orientations == 90
            metadata.StimParams.startPoint = - metadata.StimParams.right_screen;
            metadata.StimParams.endPoint = metadata.StimParams.left_screen;
        end
        metadata.StimParams.visualSpace = 'azimuth';
    case 'continousEdge_elevation'
        if metadata.StimParams.orientations == 0
            metadata.StimParams.startPoint = -metadata.StimParams.lower_screen;
            metadata.StimParams.endPoint = metadata.StimParams.upper_screen;
        elseif metadata.StimParams.orientations == 180
            metadata.StimParams.startPoint = metadata.StimParams.upper_screen;
            metadata.StimParams.endPoint =-metadata.StimParams.lower_screen;
        end
        metadata.StimParams.visualSpace = 'elevation';
    case 'Retinotopy_LR'
        Direction = regexp(metadata.StimParams.Direction,'\w*','match', 'once');
        try 
            if Direction == 'left'
                metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
                metadata.StimParams.endPoint = -metadata.StimParams.endPoint;
            end
        catch            
        end
        metadata.StimParams.visualSpace = 'azimuth';
    case 'Retinotopy_UD'
        Direction = regexp(metadata.StimParams.Direction,'\w*','match', 'once');
        try 
            if Direction == 'down'
                metadata.StimParams.startPoint = -metadata.StimParams.startPoint;
                metadata.StimParams.endPoint = -metadata.StimParams.endPoint;
            end
        catch
        end
        metadata.StimParams.visualSpace = 'elevation';
    case 'continuousEdge_left_RHField'
        metadata.StimParams.visualSpace = 'azimuth';
end

%% create stimCodes
metadata.Imaging.startTime = metadata.Imaging.time(1);
metadata.Imaging.offsetFrameTimes = metadata.Imaging.time-metadata.Imaging.startTime;
metadata.Imaging.meanCCDTime = median(diff(metadata.Imaging.offsetFrameTimes));
stimConditionOnsetTimeCorr = metadata.StimParams.stimConditionOnsetTime-metadata.Imaging.startTime;

%% make fourierMaps
metadata.Imaging.frameRate = 1/metadata.Imaging.meanCCDTime;
metadata.Imaging.stimConditionOnsetTimeInFrames = ceil(stimConditionOnsetTimeCorr*metadata.Imaging.frameRate);

%Read Imaging Data For Each Cycle
time=tic;
displayOnlineMaps = true;
harmonicsToTest = [1 2];
numberOfTrials    = length(metadata.Imaging.stimConditionOnsetTimeInFrames);
freqData   = zeros(size(tifStack,1),size(tifStack,2),length(harmonicsToTest),numberOfTrials);
Win_startFrames = floor(OnSet*metadata.Imaging.frameRate);
Win_endFrames = floor(OffSet*metadata.Imaging.frameRate);
framesPerCycle    = mode(diff(metadata.Imaging.stimConditionOnsetTimeInFrames));

for currentCycle = 1:numberOfTrials
    selectedFrames = metadata.Imaging.stimConditionOnsetTimeInFrames(currentCycle)+Win_startFrames:...
        (metadata.Imaging.stimConditionOnsetTimeInFrames(currentCycle)+framesPerCycle+Win_startFrames+Win_endFrames-1);
    temp = zeros(size(tifStack,1),size(tifStack,2),length(selectedFrames));
    if(selectedFrames(end)>size(tifStack,3)),selectedFrames = selectedFrames(1):size(tifStack,3); end
    temp(:,:,1:length(selectedFrames)) = tifStack(:,:,selectedFrames);
    temp = double(temp);
    
    % Create Temporally Encoded Maps at the Desired Frequencies using DFT
    meanTrialFrame = mean(temp,3); % Used as a Cocktail Blank
    currentRelativeHarmonic = 0;
    for currentStimFreq = harmonicsToTest
        currentRelativeHarmonic = currentRelativeHarmonic+1;
        freqData(:,:,currentRelativeHarmonic,currentCycle) = freqData(:,:,currentRelativeHarmonic,currentCycle)...
            + sum((temp-repmat(meanTrialFrame,[1 1 size(temp,3)]))...
            .* repmat(reshape(exp(2*pi*1i*currentStimFreq*([1:size(temp,3)]/framesPerCycle)),[1 1 size(temp,3)]),[size(temp,1) size(temp,2) 1]),3);
    end
    
    % Display Online Maps
    if(displayOnlineMaps && currentCycle == 1)
        h=figure; drawnow;
    end
    if(displayOnlineMaps)
        figure(h);
        imagesc(((polarMapNew(mean(freqData(:,:,1,1:currentCycle),4))))); title('1st Order Harmonic'); axis image; colormap(hsv); colorbar;
        axis equal
        set(gca,'xtick',[],'ytick',[]);
        drawnow; pause(0.25);
    end
    
    disp(['  *Finished Reading and Constructing Imaging Data for Trial ' num2str(currentCycle) ' - Time Elapsed: ' num2str(toc(time)) ' seconds']);
end
analysis.maps.polar = polarMapNew(mean(squeeze(freqData(:,:,1,:)),3),cmap,2,[0 1]);
z = mean(squeeze(freqData(:,:,1,:)),3);
analysis.maps.phaseMap = angle(z)*180/pi;
analysis.maps.angleMap = wrapTo360(analysis.maps.phaseMap)./360;
clear temp
clear tifStack
clear freqData
saveas(gcf, fullfile(saveDirectory, 'PolarMap.png'))




currentAxis=gcf;
for i = 1:currentAxis.Number
    saveas(figure(i),[saveDataDirectory 'Online map ' num2str(i) '.tif'])
end

function data = LoadImagingData(imgPath, spatialDownsamplingFactor)
    downsample = 1/spatialDownsamplingFactor;
    % Determine number of multi-page tif images are present.    
    files= dir(strcat(imgPath,'\*.tif'));
    files={files.name}';
    numberOfFiles = size(files, 1);
    
    % Determine total number of frames present, and initialize rawF array
    framesPerFile = zeros(numberOfFiles,1);
    filesToCheck = 1:numberOfFiles;
    for n = filesToCheck 
        fileName  = char(strcat(imgPath,'\',files(n)));
        imageInfo = imfinfo(fileName);
        framesPerFile(n)=numel(imageInfo);
    end
    numberOfFrames = sum(framesPerFile);
    data = zeros([imageInfo(1).Height*downsample,imageInfo(1).Width*downsample,numberOfFrames],'uint16');
    
    % Read imaging frames into MATLAB 
    frameCounter = 0;
    for n = 1:numberOfFiles
        fileName = char(strcat(imgPath,'\',files(n)));
        data(:,:,frameCounter+[1:framesPerFile(n)])=read_Tiffs(fileName,downsample,100);
        frameCounter=frameCounter+framesPerFile(n);
    end

end
function [val]= FindStimulusParam(fname,paramKey)
    fid=fopen(fname);
    try
        % get the rows for the paramKey
        fileText= textscan(fid, '%s', 'delimiter', '\n');
        fileText= fileText{1};
        fileText= strrep(fileText, ' ', ''); % delete all whitespace
        keys = strfind(fileText, strcat(paramKey, '='));
        keys= find(~cellfun(@isempty, keys));
        line = fileText(keys(1));
        line = strsplit(char(line), '#');
        line = strsplit(char(line(1)), '=');
        val= str2num(char(line(2)));
        if isempty(val)
            val = char(line(2));
        end
    catch ME
        val= '';
    end
    fclose(fid);
end
function Mask = makeMask(meanimg)
    Mean = nanmean(nanmean(meanimg,2));
    Std = std(std(meanimg));
    threshold = Mean - 0.9 * Std;
    TempMask = meanimg;
    TempMask(TempMask<threshold) = nan;
    TempMask(~isnan(TempMask)) = 1;
    TempMask(isnan(TempMask)) = 0;
    [labeledObject,nPolygons] = bwlabel(TempMask);
    MaskSize = [];
    for p = 1:nPolygons
        MaskSize(p) = sum(labeledObject(:)==p);
    end
    [~,LargestObjNr] = max(MaskSize);
    Mask = zeros(size(TempMask,1), size(TempMask,2));
    Mask(labeledObject == LargestObjNr) = 1;
end
function fig=makeFigureFullScreen(fig,constrainDimensions)
    % Makes a full-screen figure;
    if(nargin<1),fig=gcf; end
    if(nargin<2),constrainDimensions = false;  end

    screenSize = get(0,'screensize');
    set(fig,'units','normalized','outerposition',[0.05 0.05 0.9 0.9]);
    if(constrainDimensions)
        PaperPosition = get(fig,'PaperPosition');
        PaperPosition(3) = (1.25*screenSize(3)/screenSize(4))*PaperPosition(4);
        set(fig,'PaperPosition',PaperPosition);
    end
end
function HLSmap = stimulusMap(stimulusAveragedMaps, Mask)
    % Calculate max stimulus, amplitude and tuning width (circular variance of sorted map)
    NumStims = size(stimulusAveragedMaps, 3);
    dFoFmaps = zeros(size(stimulusAveragedMaps,1), size(stimulusAveragedMaps,2), size(stimulusAveragedMaps,3));
    %BSmap = mean(stimulusAveragedMaps,3);
    for stim = 1:NumStims
        I = stimulusAveragedMaps(:,:,stim);
        %dFoFmaps(:,:,stim) = (I-BSmap);
        dFoFmaps(:,:,stim) = I;
    end
    [ResponseAmplitude,PrefStim] = max(imgaussfilt(dFoFmaps), [], 3 );
    ResponseAmplitude = times(ResponseAmplitude, Mask);

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
    ResponseAmplitude_Scaled = rescale(ResponseAmplitude);
    HLSmap = zeros( size(stimulusAveragedMaps,1), size(stimulusAveragedMaps,2), 3 );            
    for y = 1:size(stimulusAveragedMaps,1)
        for x = 1:size(stimulusAveragedMaps,2)
            for c = 1:3
                if ~isnan(PrefStim(y,x)) && ResponseAmplitude(y,x) > 0
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
function rgbMap = polarMap(z)
    %calculate magnitude Map
    magnitudeMap = abs(z);
    excludingNaNs = ~isnan(magnitudeMap(:));
    highClipVal = mean(magnitudeMap(excludingNaNs(:)))+3*std(magnitudeMap(excludingNaNs(:)));
    lowClipVal  = mean(magnitudeMap(excludingNaNs(:)))-3*std(magnitudeMap(excludingNaNs(:)));
    magnitudeMap(magnitudeMap > highClipVal) = highClipVal;
    magnitudeMap(magnitudeMap < lowClipVal ) = lowClipVal;
    offsetMagnitudeMap = magnitudeMap - min(min(magnitudeMap));
    normalizedMagnitudeMap = offsetMagnitudeMap / max(max(offsetMagnitudeMap));
    
    %calculate Phase map
    phaseMap = angle(z)*180/pi;
    normalizedPhaseMap = wrapTo360(phaseMap)./360;
    normalizedPhaseMap(normalizedPhaseMap>1) = 1;
    normalizedPhaseMap(normalizedPhaseMap<0) = 0;
    
    HueMap = normalizedPhaseMap;
    SaturationMap = ones(size(z)); % normally ignored unless responseImg is specified
    ValueMap  = normalizedMagnitudeMap;
    rgbMap = hsv2rgb(cat(3,HueMap,SaturationMap,ValueMap));
end
