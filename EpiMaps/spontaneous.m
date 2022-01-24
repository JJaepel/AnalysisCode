close all
clear all

animal = 'F2537_2021-06-30 ';
expt_id =3;
sp2id = expt_id;

close all
EpiDir = 'Z:\Juliane\Data\Epi\';
Sp2Dir = 'Z:\Juliane\Data\Spike2Data\';
SaveDir = 'Z:\Juliane\Data\ImageAnalysis\';

windowStop=2;
windowStart=0;
pre=1;
field = 'rawF';
downsample = 1;


EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if expt_id > 9
%     Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id) filesep];
else
%     Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id) filesep];
end

if ~exist(saveDirectory, 'dir')
    mkdir(saveDirectory);  
end

sliceparams = struct;
sliceparams.expt_id = expt_id;
sliceparams.baseDirectory = EpiDirectory;

%% load metadata
%metadata.StimParams=Load_stimparams(Sp2dDirectory);
%metadata.Imaging=LoadFrameTimes(Sp2dDirectory);
%metadata.StimParams.path=fullfile(Sp2dDirectory);
metadata.StimParams.series=expt_id;

%% load tiffs
t0=tic;
data.rawF = readingImagingData(EpiDirectory);
toc(t0)

%% if applicaple, downsample
if downsample
    downsampleFactor = 3;
    stackSize        = size(data.rawF);
    stackSize     = floor(stackSize/downsampleFactor);
    downsampledStack = zeros(stackSize,class(data.rawF));
    for i=1:downsampleFactor
        downsampledStack = downsampledStack+data.rawF(i:downsampleFactor:downsampleFactor*stackSize(1),i:downsampleFactor:downsampleFactor*stackSize(2),i:downsampleFactor:downsampleFactor*stackSize(3))/downsampleFactor;
    end
    data.rawF = downsampledStack;
    clear downsampledStack
end

expParam.ROI =true( [size(data.rawF,1),size(data.rawF,2)]); 
expParam.rawFMeanImg = mean(data.rawF,3);
expParam.baseImg = mean(data.rawF(:,:,1:50),3);
expParam.gaussMeanImg = imgaussfilt(mean(data.rawF, 3), 4);

%% get masks for whole window as well as individual areas
[mov, mixedfilters, percent] = SVDsimple(data.rawF);
data.PCAs = mixedfilters;
figure
imagesc(data.PCAs(:,:,1))
saveas(gcf, [saveDirectory, '1stPC.png'])
figure
imagesc(data.PCAs(:,:,2))
saveas(gcf, [saveDirectory, '2dPC.png'])

mask = data.PCAs(:,:,1);
mask(mask > 0) = 0;
mask(mask < 0) = 1;
expParam.mask = logical(mask);

mask19 = data.PCAs(:,:,2);
mask19(mask19<1) = 0;
expParam.mask19 = logical(mask19);

mask19 = data.PCAs(:,:,2);
maskV1(maskV1>-200000)=0;
expParam.maskV1 = logical(maskV1);

data.filt=LowHighNormalize(double(data.rawF), expParam.mask, 1,10);
data.high = HighNormalize(doube(data.rawF), expParam.mask, 1,10);
%% get ative Frames, compute correlations of the imaging stack and show it
% Computes correlations of the imaging stack (spontaneous, response, signal, or noise). 
            % The image dimensions can be "N" dimensions, but must be organized such that they are:
            %    *Spontaneous: Should be (x,y,t). The active frames are automatically extracted. 
            %    *Response:    Collapsed from a multidimensional array (i.e. [x,y,nCond,nTrials,t])
            %                  to (x,y,n), where n is all images
            %    *Signal:      Takes in a (x,y,nCond,nTrials,t), averages along the fourth dimension,
            %    *Noise        Computes correlations along nCond.
[activeFrameStack,numberOfActiveEvents,eventOnset,eventDuration] = getActiveFrames(data.rawF,expParam);
corrTable = computeCorrelationTable(activeFrameStack,expParam.ROI);
showCorrelationStructure(corrTable,expParam, saveDirectory)