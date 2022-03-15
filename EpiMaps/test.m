%let's test out different corrleation functions and filtering

animal = 'F2537_2021-06-30';
expt_id =3;
sp2id = expt_id;


EpiDir = 'F:\Data\Epi\';
Sp2Dir = 'F:\Data\Spike2Data\';
SaveDir = 'F:\Data\ImageAnalysis\';

EpiDirectory = [EpiDir filesep animal filesep 'tseries_' num2str(expt_id) filesep];
if expt_id > 9
    Sp2dDirectory = [Sp2Dir animal filesep 't000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't000' num2str(expt_id) filesep];
else
    Sp2dDirectory = [Sp2Dir animal filesep 't0000' num2str(sp2id) filesep];
    saveDirectory = [SaveDir animal filesep 't0000' num2str(expt_id) filesep];
end

% load the raw data
load(fullfile(saveDirectory, 'rawData.mat'), 'rawF');
view_tiff_stack(rawF)

%do a quick dff based on how Alex is doing it
A=double(rawF); 
sz=size(A);
A=reshape(A,sz(1)*sz(2),sz(3));
Amean = mean(A,2); %avg at each pixel location in the image over time
A = A ./ (Amean * ones(1,sz(3))) - 1;   % F/F0 - 1 == ((F-F0)/F0);
dff=reshape(A,sz);
clear A
clear raw

%look at the normalized stack
normDa = NormalizeStack(metadata,double(rawF));

% let's remove the first PCA
[PCAfilt, mixedfilters, ~] = SVDsimple(rawF,1);

%let's look at the filtered data
ROI =true( [size(rawF,1),size(rawF,2)]);
filt=LowHighNormalizeData(dff, ROI, 1,10);
high = HighNormalizeData(dff, ROI,5);

%let's divide it into three stacks - 1)coactivity = coact, 2)
%anticorrelated = antco, 3) not active = silent
[activeFrameStack,numberOfActiveEvents,eventIndices,eventDuration] = getActiveFrames(dff);

indicesCoact = 410:490;
indicesAntCo = 1:260;
indicesSilent = 890:930;
allFrames = 1:size(dff,3);

indices = 1:(size(dff,1)*size(dff,2));

seedBasedCorr(double(normDa(:,:,eventIndices)), indices);
